# -*- coding: utf-8 -*-

"""An ORCA Frequencies (Hessian / vibrational analysis) sub-step.

Extends the Energy sub-step with ORCA's ``Freq`` (analytic, ``AnFreq``) or
``NumFreq`` (numerical) frequency calculation: the Hessian, harmonic vibrational
frequencies, IR intensities, and the thermochemistry (zero-point energy, thermal
enthalpy, entropy, and Gibbs free energy) at a chosen temperature.
"""

import logging
from pathlib import Path
import re
import textwrap

from tabulate import tabulate

import orca_step
from .energy import Energy
import seamm
from seamm_util import Q_  # noqa: F401
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
        thermochemistry."""
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
        freqs = self._parse_frequencies(text)
        if freqs is not None:
            props["frequencies"] = freqs
            n_imag = sum(1 for f in freqs if f < 0.0)
            props["n imaginary frequencies"] = n_imag
        ir = self._parse_ir_intensities(text)
        if ir is not None:
            props["IR intensities"] = ir
        props.update(self._parse_thermochemistry(text))
        props = {k: v for k, v in props.items() if v not in (None, [], {})}

        _, configuration = self.get_system_configuration(None)
        try:
            self.store_results(configuration=configuration, data=props)
        except Exception as e:  # pragma: no cover
            logger.warning(f"Could not store results: {e}")

        self._print_scalar_summary(props)
        self._report_frequencies(props)

    def _parse_frequencies(self, text):
        """The vibrational frequencies (cm^-1) from ORCA's VIBRATIONAL
        FREQUENCIES block, dropping the ~zero translations/rotations. Imaginary
        modes are reported by ORCA as negative and kept as such."""
        freqs = [
            float(v)
            for v in re.findall(r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1", text, re.M)
        ]
        vibrational = [f for f in freqs if abs(f) > 1.0]
        return vibrational or None

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
        """Zero-point energy, total enthalpy, and Gibbs free energy (all E_h)."""
        out = {}
        for key, pattern in (
            ("zero point energy", r"Zero point energy\s+\.\.\.\s+(-?\d+\.\d+)"),
            ("enthalpy", r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)"),
            ("gibbs energy", r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)"),
        ):
            m = re.findall(pattern, text)
            if m:
                out[key] = float(m[-1])
        return out

    def _report_frequencies(self, p):
        """Print the thermochemistry, note imaginary modes, and list the
        frequencies for small systems."""
        rows = []

        def add(label, value, units, fmt="{:.8f}"):
            if value is not None:
                rows.append([label, fmt.format(value), units])

        add("Zero-point energy", p.get("zero point energy"), "E_h")
        add("Total enthalpy", p.get("enthalpy"), "E_h")
        add("Gibbs free energy", p.get("gibbs energy"), "E_h")
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
