# -*- coding: utf-8 -*-

"""An ORCA geometry-optimization sub-step.

Extends the Energy sub-step by adding the ``Opt`` keyword and, after the run,
updating the configuration with the optimized geometry.
"""

import logging
from pathlib import Path
import pprint  # noqa: F401

import orca_step
from .energy import Energy
import seamm  # noqa: F401
from seamm_util.printing import FormattedText as __
import seamm_util.printing as printing

logger = logging.getLogger(__name__)
printer = printing.getPrinter("ORCA")


class Optimization(Energy):
    """A geometry optimization with ORCA.

    See Also
    --------
    TkOptimization, OptimizationParameters, Energy
    """

    def __init__(
        self, flowchart=None, title="Optimization", extension=None, logger=logger
    ):
        logger.debug(f"Creating ORCA Optimization {self}")
        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )
        self._calculation = "optimization"
        self.parameters = orca_step.OptimizationParameters()

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
        conv = P["optimization convergence"]
        text = f"Geometry optimization with ORCA at {method} ({conv})."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self, keywords=None):
        """Run the optimization by adding the Opt keyword to the energy run."""
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        opt_keywords = [P["optimization convergence"]]
        if keywords:
            opt_keywords += list(keywords)
        return super().run(keywords=opt_keywords)

    def analyze(self, indent="", P=None, data=None, **kwargs):
        """Report the energy and update the configuration with the optimized
        geometry, plus check that the optimization converged.
        """
        super().analyze(indent=indent, P=P, data=data)

        directory = Path(self.directory)

        # Convergence check from the output.
        out = directory / "orca.out"
        converged = False
        if out.exists():
            text = out.read_text()
            converged = "THE OPTIMIZATION HAS CONVERGED" in text
        if not converged:
            printer.normal(
                __(
                    "Warning: the ORCA geometry optimization did not report "
                    "convergence.",
                    indent=self.indent + 4 * " ",
                )
            )

        # Update the configuration with the optimized geometry (ORCA writes it
        # to orca.xyz).
        xyz = directory / "orca.xyz"
        if xyz.exists():
            coords = self._read_xyz_coordinates(xyz)
            _, configuration = self.get_system_configuration(None)
            if len(coords) == configuration.n_atoms:
                configuration.atoms.set_coordinates(coords, fractionals=False)
                printer.normal(
                    __(
                        "Updated the structure with the optimized geometry.",
                        indent=self.indent + 4 * " ",
                    )
                )
            else:
                logger.warning(
                    f"orca.xyz has {len(coords)} atoms, configuration has "
                    f"{configuration.n_atoms}; not updating geometry."
                )

    @staticmethod
    def _read_xyz_coordinates(path):
        """Read [[x, y, z], ...] in Angstrom from a standard .xyz file."""
        lines = path.read_text().splitlines()
        n = int(lines[0].split()[0])
        coords = []
        for line in lines[2 : 2 + n]:
            fields = line.split()
            coords.append([float(fields[1]), float(fields[2]), float(fields[3])])
        return coords
