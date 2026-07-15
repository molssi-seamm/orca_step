#!/usr/bin/env python3
"""
orca_mdi.py -- MDI engine wrapping ORCA for use as a persistent energy/force
engine behind the SEAMM MDI facility (seamm_mdi.MDIEngine) and any MDI driver.

ORCA has no in-process Python API for single-point energies and gradients, so
this engine runs the ``orca`` binary once per geometry inside a persistent
working directory. Keeping the same directory (and job name ``orca``) lets ORCA
autostart from the previous ``orca.gbw``, so each geometry's SCF starts from the
last geometry's orbitals -- the main saving for a scan of nearby structures. All
structural information (atom count, elements, coordinates) comes from the MDI
handshake; charge and multiplicity come from CLI defaults, optionally overridden
by >TOTCHARGE / >ELEC_MULT.

Note ORCA's per-geometry cost: unlike the semiempirical engines, ORCA is not
dominated by start-up, so the persistent-engine win is smaller (orbital reuse
plus not re-parsing a sub-flowchart each point). It is provided so ORCA is a
first-class MDI engine like MOPAC/xTB -- e.g. for the Dimer Builder's energy
contact search and future driver retrofits.

Usage (TCP; the driver picks the port/host):
    python orca_mdi.py \\
        -mdi "-role ENGINE -name ORCA -method TCP -port 8021 -hostname localhost" \\
        --orca /path/to/orca --method B3LYP --basis def2-SVP \\
        --charge 0 --multiplicity 1

Dependencies (all present in the SEAMM environment): pymdi, numpy, seamm_util,
and the ORCA executable (its full path passed with --orca so it can find its
sub-programs).
"""

import argparse
import logging
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from seamm_util import Q_, element_data

logger = logging.getLogger("orca-mdi")

# MDI sends coordinates in bohr; ORCA input geometries are in Angstrom. Energies
# and gradients are already in atomic units (hartree, hartree/bohr) on both
# sides, so only the geometry needs converting.
ANG_PER_BOHR = Q_(1.0, "bohr").m_as("angstrom")

# Atomic number -> element symbol, from seamm_util's element data.
_Z_TO_SYMBOL = {d["atomic number"]: sym for sym, d in element_data.items()}


# ---------------------------------------------------------------------------
# ORCA-specific helpers (pure; unit-testable without a socket or ORCA present)
# ---------------------------------------------------------------------------
def orca_input(method, basis, charge, multiplicity, symbols, coords_ang, ncores=1):
    """The ORCA input for a single-point energy + gradient (``EnGrad``).

    ``coords_ang`` is an (n, 3) array in Angstrom. A parallel run adds a
    ``%pal`` block. ORCA autostarts from an existing ``orca.gbw`` in the run
    directory, so no explicit guess keyword is needed here.
    """
    lines = [f"! {method} {basis} EnGrad"]
    if ncores and ncores > 1:
        lines.append(f"%pal nprocs {ncores} end")
    lines.append(f"* xyz {charge} {multiplicity}")
    for sym, (x, y, z) in zip(symbols, coords_ang):
        lines.append(f"{sym:2s} {x:18.10f} {y:18.10f} {z:18.10f}")
    lines.append("*")
    lines.append("")
    return "\n".join(lines)


def parse_energy(out_text):
    """The final single-point energy (hartree) from ORCA's stdout, or None."""
    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out_text)
    return float(matches[-1]) if matches else None


def orca_hessian_input(
    method, basis, charge, multiplicity, symbols, coords_ang, ncores=1
):
    """ORCA input for an analytic Hessian (``! AnFreq``), which writes an
    ``orca.hess`` file. Same geometry/charge/multiplicity handling as the
    energy+gradient input; ``AnFreq`` needs an analytic second derivative for the
    method (HF, most DFT, MP2). Methods without one should fall back to
    finite-differencing the gradient on the driver side."""
    lines = [f"! {method} {basis} AnFreq"]
    if ncores and ncores > 1:
        lines.append(f"%pal nprocs {ncores} end")
    lines.append(f"* xyz {charge} {multiplicity}")
    for sym, (x, y, z) in zip(symbols, coords_ang):
        lines.append(f"{sym:2s} {x:18.10f} {y:18.10f} {z:18.10f}")
    lines.append("*")
    lines.append("")
    return "\n".join(lines)


def parse_hessian(hess_text, natoms):
    """The Cartesian Hessian (hartree/bohr^2) as a (3n, 3n) array from an
    ``orca.hess`` file's ``$hessian`` block.

    The block is: the dimension on its own line, then column-blocks (5 columns
    each) -- a header line of column indices followed by one line per row that
    begins with the row index and lists that row's values for those columns."""
    n = 3 * natoms
    lines = hess_text.splitlines()
    start = next((i for i, ln in enumerate(lines) if ln.strip() == "$hessian"), None)
    if start is None:
        raise RuntimeError("no $hessian block in the ORCA .hess file")
    dim = int(lines[start + 1].split()[0])
    if dim != n:
        raise RuntimeError(f".hess dimension {dim} != 3*natoms ({n}).")
    H = np.zeros((n, n))
    j = start + 2
    cols_done = 0
    while cols_done < n:
        col_idx = [int(c) for c in lines[j].split()]
        j += 1
        for _ in range(n):
            parts = lines[j].split()
            row = int(parts[0])
            for k, c in enumerate(col_idx):
                H[row, c] = float(parts[1 + k])
            j += 1
        cols_done += len(col_idx)
    return H


def parse_engrad(engrad_text, natoms):
    """The Cartesian gradient (hartree/bohr) as an (natoms, 3) array from the
    contents of an ``orca.engrad`` file."""
    # The gradient block is 3*natoms floats following the header line; the file
    # interleaves '#'-comment lines, so pull every float and slice out the
    # gradient (which comes right after the single energy value).
    nums = re.findall(r"-?\d+\.\d+(?:[eE][+-]?\d+)?", engrad_text)
    # Layout: [energy, gx1, gy1, gz1, ...]; there may be trailing coordinates.
    grad = [float(x) for x in nums[1 : 1 + 3 * natoms]]
    if len(grad) != 3 * natoms:
        raise RuntimeError(
            f"orca.engrad did not contain {3 * natoms} gradient values "
            f"(found {len(grad)})."
        )
    return np.array(grad).reshape(natoms, 3)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="ORCA MDI engine (runs the orca binary)")
    p.add_argument(
        "-mdi", required=True, help="MDI initialization string passed to MDI_Init"
    )
    p.add_argument(
        "--orca",
        required=True,
        help="Full path to the orca executable (ORCA needs its full path to "
        "find its sub-programs).",
    )
    p.add_argument(
        "--method",
        required=True,
        help="The ORCA method keyword (a functional like B3LYP, or HF, MP2, "
        "DLPNO-CCSD(T), ...).",
    )
    p.add_argument(
        "--basis", default="def2-SVP", help="The orbital basis set (default def2-SVP)."
    )
    p.add_argument(
        "--auxiliary-basis",
        default="AutoAux",
        help="Auxiliary/fitting basis, appended to the method (default AutoAux). "
        "'none' omits it.",
    )
    p.add_argument("--charge", type=int, default=0, help="Total charge (default 0).")
    p.add_argument(
        "--multiplicity",
        type=int,
        default=1,
        help="Spin multiplicity 2S+1 (default 1). Overridden by >ELEC_MULT.",
    )
    p.add_argument(
        "--ncores",
        type=int,
        default=1,
        help="Cores/processes for ORCA's %%pal (default 1, serial).",
    )
    p.add_argument(
        "--scratch",
        default=None,
        help="Working directory for the ORCA runs (default: a fresh temp dir). "
        "Reused across geometries so orca.gbw persists for the SCF guess.",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s: %(message)s",
        stream=sys.stdout,
    )

    try:
        import mdi
    except ImportError:
        logger.error("mdi not found. conda install -c conda-forge pymdi")
        sys.exit(1)

    method = args.method
    aux = args.auxiliary_basis
    if aux and aux.lower() != "none":
        method = f"{method} {aux}"

    workdir = (
        Path(args.scratch)
        if args.scratch
        else Path(tempfile.mkdtemp(prefix="orca_mdi_"))
    )
    workdir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"ORCA MDI engine work dir: {workdir}")

    mdi.MDI_Init(args.mdi)
    comm = mdi.MDI_Accept_Communicator()
    mdi.MDI_Register_node("@DEFAULT")
    for _cmd in [
        "<NATOMS",
        ">NATOMS",
        "<NAME",
        ">ELEMENTS",
        ">COORDS",
        ">TOTCHARGE",
        ">ELEC_MULT",
        "SCF",
        "<ENERGY",
        "<FORCES",
        "<HESSIAN",
        "EXIT",
    ]:
        mdi.MDI_Register_command("@DEFAULT", _cmd)

    natoms = None
    atomic_numbers = None
    coords_bohr = None
    charge = args.charge
    multiplicity = args.multiplicity
    have_elements = False
    have_coords = False
    recompute = True
    energy = None
    gradient = None
    hessian = None

    def run_calculation():
        """Write the ORCA input for the current geometry, run the orca binary
        (autostarting from any orca.gbw already in the work dir), and parse the
        energy and gradient."""
        nonlocal energy, gradient
        symbols = [_Z_TO_SYMBOL[int(z)] for z in atomic_numbers]
        coords_ang = coords_bohr.reshape(natoms, 3) * ANG_PER_BOHR
        text = orca_input(
            method, args.basis, charge, multiplicity, symbols, coords_ang, args.ncores
        )
        (workdir / "orca.inp").write_text(text)
        t0 = time.perf_counter()
        with (workdir / "orca.out").open("w") as out:
            result = subprocess.run(
                [args.orca, "orca.inp"],
                cwd=workdir,
                stdout=out,
                stderr=subprocess.STDOUT,
            )
        elapsed = time.perf_counter() - t0
        if result.returncode != 0:
            raise RuntimeError(
                f"orca exited with code {result.returncode}; see {workdir}/orca.out"
            )
        out_text = (workdir / "orca.out").read_text()
        energy = parse_energy(out_text)
        if energy is None:
            raise RuntimeError(f"Could not parse an energy from {workdir}/orca.out")
        gradient = parse_engrad((workdir / "orca.engrad").read_text(), natoms)
        logger.debug(f"orca run: E = {energy:.8f} Ha in {elapsed:.2f} s")

    def run_hessian():
        """Compute the analytic Hessian for the current geometry (``! AnFreq``)
        and parse orca.hess. Separate from run_calculation because it is much
        more expensive; the result is cached until the geometry changes."""
        nonlocal hessian
        symbols = [_Z_TO_SYMBOL[int(z)] for z in atomic_numbers]
        coords_ang = coords_bohr.reshape(natoms, 3) * ANG_PER_BOHR
        text = orca_hessian_input(
            method, args.basis, charge, multiplicity, symbols, coords_ang, args.ncores
        )
        (workdir / "orca.inp").write_text(text)
        t0 = time.perf_counter()
        with (workdir / "orca.out").open("w") as out:
            result = subprocess.run(
                [args.orca, "orca.inp"],
                cwd=workdir,
                stdout=out,
                stderr=subprocess.STDOUT,
            )
        elapsed = time.perf_counter() - t0
        if result.returncode != 0:
            raise RuntimeError(
                f"orca (AnFreq) exited with code {result.returncode}; see "
                f"{workdir}/orca.out"
            )
        hessian = parse_hessian((workdir / "orca.hess").read_text(), natoms)
        logger.debug(f"orca Hessian: {3 * natoms}x{3 * natoms} in {elapsed:.2f} s")

    logger.debug("Entering MDI event loop ...")
    while True:
        command = mdi.MDI_Recv_Command(comm)
        logger.debug(f"command: {command!r}")

        if command == "<NATOMS":
            mdi.MDI_Send(natoms, 1, mdi.MDI_INT, comm)
        elif command == ">NATOMS":
            raw = mdi.MDI_Recv(1, mdi.MDI_INT, comm)
            natoms = int(np.asarray(raw).flat[0])
            atomic_numbers = np.zeros(natoms, dtype=int)
            coords_bohr = np.zeros(3 * natoms)
        elif command == "<NAME":
            mdi.MDI_Send("ORCA", mdi.MDI_NAME_LENGTH, mdi.MDI_CHAR, comm)
        elif command == ">ELEMENTS":
            raw = mdi.MDI_Recv(natoms, mdi.MDI_INT, comm)
            atomic_numbers[:] = np.asarray(raw)
            have_elements = True
            recompute = True
            hessian = None
        elif command == ">COORDS":
            raw = mdi.MDI_Recv(3 * natoms, mdi.MDI_DOUBLE, comm)
            coords_bohr[:] = np.asarray(raw)
            have_coords = True
            recompute = True
            hessian = None
        elif command == ">TOTCHARGE":
            raw = mdi.MDI_Recv(1, mdi.MDI_DOUBLE, comm)
            charge = int(round(float(np.asarray(raw).flat[0])))
            recompute = True
            hessian = None
        elif command == ">ELEC_MULT":
            raw = mdi.MDI_Recv(1, mdi.MDI_INT, comm)
            multiplicity = int(np.asarray(raw).flat[0])
            recompute = True
            hessian = None
        elif command in ("SCF", "<ENERGY", "<FORCES"):
            if recompute:
                if not (have_elements and have_coords):
                    raise RuntimeError(
                        f"{command} requested before >ELEMENTS and >COORDS were "
                        "received."
                    )
                run_calculation()
                recompute = False
            if command == "<ENERGY":
                mdi.MDI_Send(energy, 1, mdi.MDI_DOUBLE, comm)
            elif command == "<FORCES":
                forces = (-gradient).ravel()  # force = -dE/dx, in hartree/bohr
                mdi.MDI_Send(forces, 3 * natoms, mdi.MDI_DOUBLE, comm)
        elif command == "<HESSIAN":
            # The analytic Cartesian Hessian (hartree/bohr^2), 3N x 3N row-major.
            # Custom command: a driver introspects @COMMANDS and uses this when
            # present, else finite-differences <FORCES itself.
            if not (have_elements and have_coords):
                raise RuntimeError(
                    "<HESSIAN requested before >ELEMENTS and >COORDS were received."
                )
            if hessian is None:
                run_hessian()
            mdi.MDI_Send(
                np.asarray(hessian).ravel(),
                9 * natoms * natoms,
                mdi.MDI_DOUBLE,
                comm,
            )
        elif command == "EXIT":
            logger.debug("EXIT -- shutting down.")
            break
        else:
            logger.warning(f"unrecognised command {command!r}")

    logger.debug("Done.")


if __name__ == "__main__":
    main()
