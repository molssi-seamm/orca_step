# -*- coding: utf-8 -*-

"""An ORCA Frequencies (Hessian / vibrational analysis) sub-step.

Extends the Energy sub-step with ORCA's ``Freq`` (analytic, ``AnFreq``) or
``NumFreq`` (numerical) frequency calculation: the Hessian, harmonic vibrational
frequencies, IR intensities, and the thermochemistry (zero-point energy, thermal
enthalpy, entropy, and Gibbs free energy) at a chosen temperature.
"""

import csv
import logging
import math
from pathlib import Path
import re
import textwrap

import numpy as np
from tabulate import tabulate

import orca_step
from .energy import Energy
import seamm
from seamm_util import Q_
from seamm_util.printing import FormattedText as __
import seamm_util.printing as printing

logger = logging.getLogger(__name__)
printer = printing.getPrinter("ORCA")

# Conversion from a mass-weighted Cartesian force constant, in
# Hartree/(bohr**2 * u), to a harmonic wavenumber in cm**-1:
#     nu = sqrt(lambda) * _CM_PER_SQRT_AU
# Derived from CODATA constants as sqrt(E_h / (a0**2 * u)) / (2*pi*c) ~ 5140.49.
_E_H = 4.3597447222071e-18  # Hartree, J
_A0 = 5.29177210903e-11  # Bohr radius, m
_AMU = 1.66053906660e-27  # atomic mass unit, kg
_C_CM = 2.99792458e10  # speed of light, cm/s
_CM_PER_SQRT_AU = math.sqrt(_E_H / (_A0**2 * _AMU)) / (2 * math.pi * _C_CM)


class Frequencies(Energy):
    """A vibrational-frequency (Hessian) calculation with ORCA.

    See Also
    --------
    TkFrequencies, FrequenciesParameters, Energy
    """

    def __init__(
        self, flowchart=None, title="Frequencies", extension=None, logger=logger
    ):
        logger.debug(f"Creating ORCA Frequencies {self}")
        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )
        self._calculation = "frequencies"
        self.parameters = orca_step.FrequenciesParameters()

    def _extrapolating(self, P):
        """Never extrapolate: a CBS-extrapolated energy has no Hessian."""
        return False

    def description_text(self, P=None):
        if not P:
            P = self.parameters.values_to_dict()
        how = "numerical" if P.get("second derivatives") == "numerical" else "analytic"
        text = (
            f"Vibrational frequencies with ORCA at {self._level_of_theory_text(P)} "
            f"({how} Hessian), with thermochemistry."
        )
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self, keywords=None):
        """Add ORCA's Freq/NumFreq keyword to the energy run."""
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        freq_kw = "NumFreq" if P.get("second derivatives") == "numerical" else "Freq"
        kws = [freq_kw]
        if keywords:
            kws += list(keywords)
        return super().run(keywords=kws)

    def extra_input(self, P):
        """Add the ``%freq Temp`` block (thermochemistry temperature) to the
        Energy sub-step's extra ORCA input."""
        blocks, files = super().extra_input(P)
        temp = P.get("temperature", 298.15)
        temperature = temp.m_as("K") if hasattr(temp, "m_as") else float(temp)
        freq_block = f"%freq Temp {temperature:.4f} end"
        blocks = f"{blocks}\n{freq_block}" if blocks.strip() else freq_block
        return blocks, files

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    def analyze(self, indent="", P=None, data=None, **kwargs):
        """Store and report the energy plus the frequencies, IR intensities, and
        thermochemistry, honoring the standard structure-handling options."""
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )
        directory = Path(self.directory)
        props = self._parse_properties(directory)
        if "energy" not in props and data and data.get("energy") is not None:
            props["energy"] = data["energy"]

        text = (
            (directory / "orca.out").read_text()
            if (directory / "orca.out").exists()
            else ""
        )

        _, initial_configuration = self.get_system_configuration(None)

        all_freqs = self._parse_frequencies(text)
        ir = self._parse_ir_intensities(text)
        if all_freqs is not None:
            vibrational, _ = self._classify_frequencies(
                all_freqs, initial_configuration
            )
            props["frequencies"] = vibrational
            props["n imaginary frequencies"] = sum(1 for f in vibrational if f < 0.0)
            # ORCA projects the translations/rotations to exactly 0.00 in its
            # printed frequencies, so the genuine numerical residual (a gauge of
            # the Hessian's accuracy) is obtained by diagonalizing the raw,
            # un-projected mass-weighted Hessian from orca.hess.
            residual = self._zero_mode_residual(directory, initial_configuration)
            if residual is not None:
                props["largest zero-mode frequency"] = residual
            if ir is not None:
                props["IR intensities"] = ir
            # Write the frequencies (and IR intensities) to a CSV for easy access.
            self._write_frequencies_csv(directory, vibrational, ir)
            # Write a graph of the IR spectrum (sticks + a broadened trace).
            if ir is not None and len(ir) == len(vibrational):
                try:
                    self._plot_ir_spectrum(directory, vibrational, ir)
                except Exception as e:  # pragma: no cover
                    logger.warning(f"Could not create the IR-spectrum graph: {e}")

        props.update(self._parse_thermochemistry(text))
        props = {k: v for k, v in props.items() if v not in (None, [], {})}

        # Store the results/properties per the structure-handling options. A
        # frequency calculation does not change the geometry, so a new
        # configuration simply carries a copy of the current structure.
        handling = P.get("structure handling", "Overwrite the current configuration")
        if handling == "Discard the structure":
            system, configuration = None, initial_configuration
        else:
            system, configuration = self.get_system_configuration(
                P, same_as=initial_configuration
            )
        try:
            self.store_results(configuration=configuration, data=props)
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not store results: {e}")

        self._print_scalar_summary(props)
        self._report_frequencies(props)
        if handling != "Discard the structure":
            printer.normal(
                __(
                    seamm.standard_parameters.set_names(system, configuration, P),
                    indent=self.indent + 4 * " ",
                )
            )

    def _parse_frequencies(self, text):
        """All 3N harmonic frequencies (cm^-1) from ORCA's VIBRATIONAL
        FREQUENCIES block, in ascending order (including the ~zero
        translations/rotations). Imaginary modes are reported by ORCA as
        negative and kept as such."""
        freqs = [
            float(v)
            for v in re.findall(r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1", text, re.M)
        ]
        return freqs or None

    def _classify_frequencies(self, all_freqs, configuration):
        """Split ORCA's 3N frequencies into the true vibrational modes and the 5
        (linear) or 6 (non-linear) nominally-zero translations/rotations.

        Returns ``(vibrational, max_zero)`` where ``vibrational`` are the
        remaining modes (ascending, imaginary modes kept negative) and
        ``max_zero`` is the largest ``|frequency|`` among the trans/rot modes --
        a gauge of the numerical accuracy of the Hessian (it should be near
        zero).
        """
        n_tr = 5 if self._is_linear(configuration) else 6
        n_tr = min(n_tr, len(all_freqs))
        if n_tr == 0:
            return list(all_freqs), None
        # The translations/rotations are the n_tr modes closest to zero, so a
        # genuine imaginary mode (a transition state) stays with the vibrations.
        order = sorted(range(len(all_freqs)), key=lambda i: abs(all_freqs[i]))
        tr = set(order[:n_tr])
        vibrational = [f for i, f in enumerate(all_freqs) if i not in tr]
        max_zero = max(abs(all_freqs[i]) for i in tr)
        return vibrational, max_zero

    @staticmethod
    def _is_linear(configuration):
        """Whether the molecule is linear (so it has 5 rather than 6 zero
        modes), from the rank of the centered coordinates."""
        pts = np.asarray(
            configuration.atoms.get_coordinates(fractionals=False), dtype=float
        )
        if len(pts) <= 2:
            return True
        s = np.linalg.svd(pts - pts.mean(axis=0), compute_uv=False)
        return bool(s[1] < 1.0e-3 * s[0])

    def _zero_mode_residual(self, directory, configuration):
        """The largest ``|frequency|`` (cm**-1) among the 5 (linear) or 6
        (non-linear) translational/rotational modes, from the **raw,
        un-projected** mass-weighted Hessian in ``orca.hess``.

        ORCA projects these modes to exactly ``0.00`` in its printed
        frequencies, so this diagonalization is the only way to see the genuine
        numerical residual -- a useful gauge of the Hessian's accuracy. Returns
        ``None`` when ``orca.hess`` is not available.
        """
        parsed = self._parse_hess_file(Path(directory) / "orca.hess")
        if parsed is None:
            return None
        hessian, masses = parsed
        # Mass-weight the Cartesian Hessian: divide each element by
        # sqrt(m_i * m_j) with the mass repeated over the three coordinates.
        m = np.repeat(masses, 3)
        if hessian.shape != (m.size, m.size):
            return None
        mass_weighted = hessian / np.sqrt(np.outer(m, m))
        evals = np.linalg.eigvalsh(mass_weighted)
        # A negative eigenvalue is an imaginary mode; keep the sign so the
        # magnitude (not the sign) drives the "closest to zero" selection.
        freqs = np.sign(evals) * np.sqrt(np.abs(evals)) * _CM_PER_SQRT_AU
        n_tr = 5 if self._is_linear(configuration) else 6
        n_tr = min(n_tr, len(freqs))
        if n_tr == 0:
            return None
        residual = sorted(freqs, key=abs)[:n_tr]
        return max(abs(float(f)) for f in residual)

    @staticmethod
    def _parse_hess_file(path):
        """Read ORCA's ``.hess`` file: the Cartesian Hessian (Hartree/bohr**2,
        as a ``(3N, 3N)`` array) and the atomic masses (u). Returns
        ``(hessian, masses)`` or ``None`` if the file or a block is missing."""
        if not Path(path).exists():
            return None
        lines = Path(path).read_text().splitlines()

        def find(tag):
            for k, line in enumerate(lines):
                if line.strip() == tag:
                    return k
            return None

        # Masses from the $atoms block: "<symbol> <mass> <x> <y> <z>".
        i = find("$atoms")
        j = find("$hessian")
        if i is None or j is None:
            return None
        n_atoms = int(lines[i + 1].split()[0])
        masses = np.array([float(lines[i + 2 + k].split()[1]) for k in range(n_atoms)])

        # The Hessian is written in blocks of (up to) 5 columns: a header line
        # of column indices, then one line per row ("<row> v0 v1 ... v4").
        dim = int(lines[j + 1].split()[0])
        hessian = np.zeros((dim, dim))
        p = j + 2
        while p < len(lines):
            header = lines[p].split()
            if not header or not all(x.lstrip("-").isdigit() for x in header):
                break
            cols = [int(x) for x in header]
            p += 1
            for _ in range(dim):
                row = lines[p].split()
                r = int(row[0])
                for c, value in zip(cols, row[1:]):
                    hessian[r][c] = float(value)
                p += 1
        return hessian, masses

    @staticmethod
    def _write_frequencies_csv(directory, frequencies, intensities):
        """Write the vibrational frequencies (and IR intensities, if available)
        to ``frequencies.csv`` in the step directory."""
        have_ir = intensities is not None and len(intensities) == len(frequencies)
        with open(Path(directory) / "frequencies.csv", "w", newline="") as fd:
            w = csv.writer(fd)
            header = ["Mode", "Frequency (cm^-1)"]
            if have_ir:
                header.append("IR intensity (km/mol)")
            w.writerow(header)
            for i, f in enumerate(frequencies, start=1):
                row = [i, f"{f:.2f}"]
                if have_ir:
                    row.append(f"{intensities[i - 1]:.2f}")
                w.writerow(row)

    def _plot_ir_spectrum(self, directory, frequencies, intensities, fwhm=15.0):
        """Write ``IR_spectrum.graph`` with the IR spectrum as a stick trace plus
        a Lorentzian-broadened trace (FWHM ``fwhm`` cm^-1) that mimics an
        experimental spectrum. Only the real (positive) modes are plotted."""
        modes = [
            (f, a)
            for f, a in zip(frequencies, intensities)
            if f > 0.0 and a is not None
        ]
        if not modes:
            return
        freqs = np.array([f for f, _ in modes], dtype=float)
        ints = np.array([a for _, a in modes], dtype=float)

        # A wavenumber grid spanning the modes with a little padding.
        lo = max(0.0, float(freqs.min()) - 200.0)
        hi = float(freqs.max()) + 200.0
        grid = np.arange(lo, hi + 1.0, 1.0)

        # Sum of Lorentzians, each normalized so its peak height is the mode's
        # intensity (so the broadened curve is on the same scale as the sticks).
        half = fwhm / 2.0
        broadened = np.zeros_like(grid)
        for f, a in zip(freqs, ints):
            broadened += a * half**2 / ((grid - f) ** 2 + half**2)

        figure = self.create_figure(
            module_path=(self.__module__.split(".")[0], "seamm"),
            template="line.graph_template",
            title="IR spectrum",
        )
        plot = figure.add_plot("IR")
        x_axis = plot.add_axis("x", label="Wavenumber (cm<sup>-1</sup>)")
        y_axis = plot.add_axis("y", label="IR intensity (km/mol)", anchor=x_axis)
        x_axis.anchor = y_axis

        # The broadened trace first so the sticks sit on top.
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="broadened",
            x=grid.tolist(),
            xlabel="wavenumber",
            xunits="1/cm",
            y=broadened.tolist(),
            ylabel="intensity",
            yunits="km/mol",
            color="#4dbd74",
        )

        # The stick spectrum: a vertical segment (0 -> intensity) per mode,
        # drawn as one line trace broken by None between sticks.
        xs = []
        ys = []
        for f, a in zip(freqs, ints):
            xs += [float(f), float(f), None]
            ys += [0.0, float(a), None]
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="sticks",
            x=xs,
            xlabel="wavenumber",
            xunits="1/cm",
            y=ys,
            ylabel="intensity",
            yunits="km/mol",
            color="black",
        )

        figure.grid_plots("IR")
        figure.write_file(Path(directory) / "IR_spectrum.graph")

        # Extra formats (png, etc.) if the parent step's options request them.
        options = getattr(getattr(self, "parent", None), "options", None) or {}
        formats = options.get("graph_formats")
        if formats:
            if isinstance(formats, str):
                import shlex

                formats = shlex.split(formats)
            for _format in formats:
                figure.write_file(
                    Path(directory) / f"IR_spectrum.{_format}",
                    width=int(options.get("graph_width", 1024)),
                    height=int(options.get("graph_height", 1024)),
                )

    def _parse_ir_intensities(self, text):
        """The IR intensities (km/mol) from the IR SPECTRUM block, in mode order
        (aligned with the vibrational frequencies)."""
        m = self._last(r"IR SPECTRUM\s*\n-+\n(.*?)\n\s*\n", text, re.DOTALL)
        if not m:
            return None
        intensities = []
        for line in m.group(1).splitlines():
            # Mode line: "  6:   1790.72   0.015824   79.97  ..."; Int is column 3.
            mm = re.match(r"\s*\d+:\s+-?\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)", line)
            if mm:
                intensities.append(float(mm.group(1)))
        return intensities or None

    def _parse_thermochemistry(self, text):
        """Zero-point energy, total enthalpy, and Gibbs free energy, converted
        from ORCA's E_h to kJ/mol (SEAMM's SI-based default energy unit)."""
        out = {}
        for key, pattern in (
            ("zero point energy", r"Zero point energy\s+\.\.\.\s+(-?\d+\.\d+)"),
            ("enthalpy", r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)"),
            ("gibbs energy", r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)"),
        ):
            m = re.findall(pattern, text)
            if m:
                out[key] = Q_(float(m[-1]), "E_h").m_as("kJ/mol")
        return out

    def _report_frequencies(self, p):
        """Print the thermochemistry, note imaginary modes, and list the
        frequencies for small systems."""
        rows = []

        def add(label, value, units, fmt="{:.3f}"):
            if value is not None:
                rows.append([label, fmt.format(value), units])

        add("Zero-point energy", p.get("zero point energy"), "kJ/mol")
        add("Total enthalpy", p.get("enthalpy"), "kJ/mol")
        add("Gibbs free energy", p.get("gibbs energy"), "kJ/mol")
        if rows:
            tmp = tabulate(
                rows,
                headers=["Thermochemistry", "Value", "Units"],
                tablefmt="rounded_outline",
                colalign=("left", "right", "left"),
                disable_numparse=True,
            )
            printer.normal("")
            printer.normal(textwrap.indent(tmp, self.indent + 7 * " "))

        freqs = p.get("frequencies")
        if freqs:
            n_imag = p.get("n imaginary frequencies", 0)
            printer.normal(
                __(
                    f"Found {len(freqs)} vibrational modes"
                    + (
                        f", including {n_imag} imaginary "
                        f"(negative) frequenc{'y' if n_imag == 1 else 'ies'} -- "
                        "the structure is not a minimum."
                        if n_imag
                        else " (all real; the structure is a minimum)."
                    ),
                    indent=self.indent + 4 * " ",
                )
            )

            # A table of the frequencies (and IR intensities, if available).
            ir = p.get("IR intensities")
            have_ir = ir is not None and len(ir) == len(freqs)
            headers = ["Mode", "Frequency (cm^-1)"]
            colalign = ["right", "right"]
            if have_ir:
                headers.append("IR intensity (km/mol)")
                colalign.append("right")
            table = []
            for i, f in enumerate(freqs, start=1):
                row = [i, f"{f:.2f}"]
                if have_ir:
                    row.append(f"{ir[i - 1]:.2f}")
                table.append(row)
            tmp = tabulate(
                table,
                headers=headers,
                tablefmt="rounded_outline",
                colalign=colalign,
                disable_numparse=True,
            )
            printer.normal("")
            printer.normal(textwrap.indent(tmp, self.indent + 7 * " "))

        max_zero = p.get("largest zero-mode frequency")
        if max_zero is not None:
            printer.normal(
                __(
                    "The largest of the nominally-zero translation/rotation "
                    f"frequencies is {max_zero:.2f} cm**-1 (from the raw, "
                    "un-projected Hessian; ORCA zeroes these in its printed "
                    "frequencies). It should be small -- it gauges the numerical "
                    "accuracy of the Hessian.",
                    indent=self.indent + 4 * " ",
                )
            )
