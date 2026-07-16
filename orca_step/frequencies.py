# -*- coding: utf-8 -*-

"""An ORCA Frequencies (Hessian / vibrational analysis) sub-step.

Extends the Energy sub-step with ORCA's ``Freq`` (analytic, ``AnFreq``) or
``NumFreq`` (numerical) frequency calculation: the Hessian, harmonic vibrational
frequencies, IR intensities, and the thermochemistry (zero-point energy, thermal
enthalpy, entropy, and Gibbs free energy) at a chosen temperature.
"""

import csv
import logging
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
        use_mc = P["use model chemistry"]
        if not isinstance(use_mc, bool):
            use_mc = use_mc == "yes"
        if use_mc:
            method = "the model chemistry"
        else:
            m, basis = self._resolve_method_basis(P)
            method = f"{m}/{basis}"
        how = "numerical" if P.get("second derivatives") == "numerical" else "analytic"
        text = (
            f"Vibrational frequencies with ORCA at {method} ({how} Hessian), with "
            "thermochemistry."
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
            vibrational, max_zero = self._classify_frequencies(
                all_freqs, initial_configuration
            )
            props["frequencies"] = vibrational
            props["n imaginary frequencies"] = sum(1 for f in vibrational if f < 0.0)
            if max_zero is not None:
                props["largest zero-mode frequency"] = max_zero
            if ir is not None:
                props["IR intensities"] = ir
            # Write the frequencies (and IR intensities) to a CSV for easy access.
            self._write_frequencies_csv(directory, vibrational, ir)

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

        max_zero = p.get("largest zero-mode frequency")
        if max_zero is not None:
            printer.normal(
                __(
                    "The largest of the nominally-zero translation/rotation "
                    f"frequencies is {max_zero:.2f} cm**-1, a gauge of the "
                    "numerical accuracy of the Hessian (it should be near zero).",
                    indent=self.indent + 4 * " ",
                )
            )
