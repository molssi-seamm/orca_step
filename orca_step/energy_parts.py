# -*- coding: utf-8 -*-

"""The parts of a correlated energy from ORCA: the reference (HF) energy, the CABS
singles of F12, the correlation energies and the triples, for composite recipes
that combine parts of calculations (composite_step) and for the fragment store,
which keeps them with the energy (its ``byproducts``).

The layout (agreed with the thermochemistry science, 2026-11-12)::

    {"parts version": 1, "units": "E_h",
     "parts": {name: value},              # energies, E_h, as ORCA gives them
     "total": [names],                    # the parts that sum to the energy
     "factors": {"(T) scale f": value},   # dimensionless
     "reference": "RHF" | "UHF" | "QRO",
     "frozen core orbitals": n, "relativistic": "none" | "DKH2" | "X2C",
     "diagnostics": {"T1": ..., "<S**2>": ...},
     "sources": {name: "property file" | "output" | "computed"}}

The parts of each kind of run:

=====================  ===============================================================
MP2                    HF, MP2 correlation
RI-MP2                 HF, RI-MP2 correlation
CCSD(T)                HF, CCSD correlation, (T); MP2 correlation (closed-shell RHF)
CCSD(T)-F12D/RI        HF, CABS singles, CCSD correlation, F12 correction, (T); MP2
                       correlation (closed-shell RHF); factor (T) scale f
DLPNO-CCSD(T)          HF, DLPNO-MP2 correlation, DLPNO-CCSD correlation, DLPNO-(T0)
DLPNO-CCSD(T)-F12D     HF, CABS singles, DLPNO-CCSD correlation, F12 correction,
                       DLPNO-(T0), DLPNO-MP2 correlation; factor (T) scale f
                       (computed)
=====================  ===============================================================

"HF" is the energy of the reference determinant the correlation is built on: the
SCF energy, or for quasi-restricted orbitals (QRO) the QRO determinant's energy,
not the UHF one. For F12 it excludes the CABS singles. "(T)" is never scaled; f =
(CCSD correlation + F12 correction) / CCSD correlation, as ORCA prints it for
canonical F12 and computed the same way for DLPNO-F12. The DLPNO correlation is
everything ORCA's total includes (weak pairs, PNO corrections). The MP2 energy a
CC run prints is reported only for closed-shell RHF; for UHF or QRO it is not a
standalone method and is left out.

The parts must reproduce ORCA's final single-point energy within 2e-9 E_h (the
splits are printed to 1e-9 E_h); otherwise no parts are given.
"""

import logging
import re

logger = logging.getLogger(__name__)

PARTS_VERSION = 1
#: The parts must sum to the final energy within this, E_h
SUM_TOLERANCE = 2e-9

_NUMBER = r"(-?\d+\.\d+(?:[eE][-+]?\d+)?)"


def property_blocks(text):
    """{block name: {key: value}} from an ORCA ``.property.txt``: numbers (the
    last number of an array of one value), integers and strings."""
    blocks = {}
    current = None
    key = None
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("$"):
            name = stripped[1:].strip()
            if name == "End":
                current = None
            else:
                current = blocks.setdefault(name, {})
            key = None
            continue
        if current is None:
            continue
        if stripped.startswith("&"):
            m = re.match(r"&(\S+)\s+\[&Type \"(\w+)\"[^\]]*\]\s*(.*)$", stripped)
            if not m:
                continue
            key, kind, rest = m.groups()
            if kind in ("Double", "Integer", "Boolean"):
                value = rest.split()[0] if rest else None
                if kind == "Double":
                    current[key] = float(value)
                elif kind == "Integer":
                    current[key] = int(value)
                else:
                    current[key] = value == "true"
                key = None
            elif kind == "String":
                current[key] = rest.strip().strip('"')
                key = None
            continue
        if key is not None:
            m = re.match(r"^0\s+" + _NUMBER + r"\s*$", stripped)
            if m:
                current[key] = float(m.group(1))
                key = None
    return blocks


def _last(pattern, text):
    found = re.findall(pattern, text, flags=re.MULTILINE)
    return float(found[-1]) if found else None


def energy_parts(output, properties_text):
    """The parts of the run's energy, or None (with the reason logged) when the
    run has none or they cannot be given safely.

    Parameters
    ----------
    output : str
        ``orca.out``.
    properties_text : str
        ``orca.property.txt``.
    """
    if not output or not properties_text:
        return None
    blocks = property_blocks(properties_text)
    final = blocks.get("Single_Point_Data", {}).get("FinalEnergy")
    scf = blocks.get("SCF_Energy", {}).get("totalEnergy")
    if final is None or scf is None:
        return None

    hf_type = re.search(r"Hartree-Fock type\s+HFTyp\s+\.+\s+(\w+)", output)
    hf_type = hf_type.group(1).upper() if hf_type else "RHF"
    if hf_type not in ("RHF", "UHF"):
        logger.info(f"No energy parts for a {hf_type} reference")
        return None
    reference = hf_type
    if re.search(r"UseQROs\s+true", output, flags=re.IGNORECASE):
        reference = "QRO"
    upper_input = "\n".join(
        line for line in output.splitlines() if re.match(r"\|\s*\d+>", line)
    ).upper()
    relativistic = "none"
    if re.search(r"\bDKH2?\b", upper_input):
        relativistic = "DKH2"
    elif re.search(r"\bX2C\b", upper_input):
        relativistic = "X2C"

    parts = {}
    sources = {}
    factors = {}
    diagnostics = {}
    frozen = None

    mdci = blocks.get("MDCI_Energies")
    mp2 = blocks.get("MP2_Energies")
    if mdci:
        f12 = "F12-HF basis set correction" in output
        dlpno = "DLPNO SETTINGS" in output or "E(L-MP2)" in output
        ref = mdci.get("refEnergy")
        corr = mdci.get("corrEnergy")
        triples = mdci.get("triplesEnergy")
        if None in (ref, corr):
            return None
        if mdci.get("numOfEl") is not None and mdci.get("numOfCorrEl") is not None:
            frozen = (mdci["numOfEl"] - mdci["numOfCorrEl"]) // 2
        prefix = "DLPNO-" if dlpno else ""
        t_name = None
        if triples is not None and triples != 0.0:
            if dlpno:
                t_name = "DLPNO-(T1)" if "(T1)" in upper_input else "DLPNO-(T0)"
            else:
                t_name = "(T)"
        if f12:
            cabs = _last(r"^F12-HF basis set correction\s+\.+\s+" + _NUMBER, output)
            correction = _last(
                r"^F12-correction to the correlation energy\s+\.+\s+" + _NUMBER,
                output,
            )
            if None in (cabs, correction):
                return None
            # The CCSD correlation is all of the total correlation but the F12
            # correction and the triples, so that it includes whatever ORCA's total
            # does (for DLPNO, more than its printed 'uncorrected' value)
            uncorrected = corr - (triples or 0.0) - correction
            # MDCI's reference energy is HF + CABS singles for F12
            parts["HF"] = ref - cabs
            sources["HF"] = "property file"
            parts["CABS singles"] = cabs
            sources["CABS singles"] = "output"
            parts[f"{prefix}CCSD correlation"] = uncorrected
            sources[f"{prefix}CCSD correlation"] = "computed"
            parts["F12 correction"] = correction
            sources["F12 correction"] = "output"
            printed = _last(
                r"^f = CCSD \(with F12\)/ CCSD \(without F12\)\s+\.+\s+" + _NUMBER,
                output,
            )
            if printed is not None:
                factors["(T) scale f"] = printed
                sources["(T) scale f"] = "output"
            else:
                factors["(T) scale f"] = (uncorrected + correction) / uncorrected
                sources["(T) scale f"] = "computed"
        else:
            parts["HF"] = ref
            sources["HF"] = "property file"
            ccsd = corr - (triples or 0.0)
            parts[f"{prefix}CCSD correlation"] = ccsd
            sources[f"{prefix}CCSD correlation"] = "property file"
        if t_name is not None:
            parts[t_name] = triples
            sources[t_name] = "property file"
        # The MP2 energy MDCI prints first
        if dlpno:
            value = _last(
                r"^E\(L-MP2\) including corrections\s+\.+\s+" + _NUMBER, output
            )
            if value is not None:
                parts["DLPNO-MP2 correlation"] = value
                sources["DLPNO-MP2 correlation"] = "output"
        elif reference == "RHF":
            value = _last(r"^E\(MP2\)\s+\.+\s+" + _NUMBER, output)
            if value is not None:
                parts["MP2 correlation"] = value
                sources["MP2 correlation"] = "output"
        t1 = _last(r"^T1 diagnostic\s+\.+\s+" + _NUMBER, output)
        if t1 is not None:
            diagnostics["T1"] = t1
        # The sum excludes the MP2 by-product, which is not part of the total
        summed = [k for k in parts if "MP2" not in k]
    elif mp2:
        ref = mp2.get("refEnergy")
        corr = mp2.get("corrEnergy")
        if None in (ref, corr):
            return None
        name = (
            "RI-MP2 correlation"
            if "RI-MP2 CORRELATION ENERGY" in output
            else ("MP2 correlation")
        )
        parts["HF"] = ref
        sources["HF"] = "property file"
        parts[name] = corr
        sources[name] = "property file"
        ncore = re.search(r"Freezing NCore=(\d+)", output)
        if ncore:
            frozen = int(ncore.group(1)) // 2
        summed = list(parts)
    else:
        return None

    s2 = _last(r"Expectation value of <S\*\*2>\s*:\s*" + _NUMBER, output)
    if s2 is not None:
        diagnostics["<S**2>"] = s2

    total = sum(parts[k] for k in summed)
    if abs(total - final) > SUM_TOLERANCE:
        logger.warning(
            f"The energy parts sum to {total:.10f} E_h, not the final energy "
            f"{final:.10f}; no parts are given."
        )
        return None

    return {
        "parts version": PARTS_VERSION,
        "units": "E_h",
        "parts": parts,
        "total": summed,
        "factors": factors,
        "reference": reference,
        "frozen core orbitals": frozen,
        "relativistic": relativistic,
        "diagnostics": diagnostics,
        "sources": sources,
    }
