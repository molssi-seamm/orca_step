# -*- coding: utf-8 -*-

"""An ORCA single-point energy sub-step."""

import logging
import pprint  # noqa: F401

import orca_step
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("ORCA")


class Energy(orca_step.ORCABase):
    """A single-point energy with ORCA.

    See Also
    --------
    TkEnergy, EnergyParameters, Optimization
    """

    def __init__(self, flowchart=None, title="Energy", extension=None, logger=logger):
        """Initialize the Energy sub-step."""
        logger.debug(f"Creating ORCA Energy {self}")

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._calculation = "energy"
        self._metadata = orca_step.metadata
        self.parameters = orca_step.EnergyParameters()

    def description_text(self, P=None):
        """Describe what this sub-step will do."""
        if not P:
            P = self.parameters.values_to_dict()

        use_mc = P["use model chemistry"]
        if not isinstance(use_mc, bool):
            use_mc = use_mc == "yes"

        if use_mc:
            text = (
                "Single-point energy with ORCA using the model chemistry from a "
                "preceding Model Chemistry step."
            )
        else:
            text = f"Single-point energy with ORCA at {P['method']}/{P['basis']}."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def keyword_line(self, P):
        """Build the ORCA '!' keyword line.

        If ``use model chemistry`` is set, the method (and basis, if given) come
        from the global ``_model_chemistry`` variable defined by a preceding
        Model Chemistry step; ORCA raises a clear error if it cannot provide that
        model chemistry. Otherwise the explicit method/basis parameters are used.
        The auxiliary basis and extra keywords always come from this node.
        """
        use_mc = P["use model chemistry"]
        if not isinstance(use_mc, bool):
            use_mc = use_mc == "yes"

        if use_mc:
            method, basis = self._method_basis_from_model_chemistry(P)
        else:
            method, basis = P["method"], P["basis"]

        keywords = [method, basis]
        aux = P["auxiliary basis"]
        if aux and aux.lower() != "none":
            keywords.append(aux)
        extra = P["extra keywords"].strip()
        if extra:
            keywords.append(extra)
        return " ".join(k for k in keywords if k)

    def _method_basis_from_model_chemistry(self, P):
        """Resolve (method, basis) from the global model chemistry, erroring
        clearly if ORCA cannot provide it.
        """
        if not self.variable_exists("_model_chemistry"):
            raise RuntimeError(
                "'Use the global model chemistry' is selected, but no Model "
                "Chemistry step has defined one (the '_model_chemistry' variable "
                "is not set). Add a Model Chemistry step before this ORCA step, "
                "or turn off 'use model chemistry' and set the method and basis "
                "explicitly."
            )
        mc = self.get_variable("_model_chemistry")
        level = mc.get("level", "") or mc.get("method", "")
        owner = (mc.get("owner") or "").strip()
        mtype = (mc.get("type") or "").strip()
        method = (mc.get("method") or "").strip()
        basis = mc.get("basis") or ""

        # ORCA can evaluate an ORCA-owned or program-agnostic level of theory,
        # but not another program's (e.g. MOPAC/xTB/LAMMPS) or a non-QC type.
        if owner not in ("", "ORCA"):
            raise RuntimeError(
                f"The model chemistry '{level}' is owned by '{owner}', which "
                "ORCA cannot provide. Choose an ORCA (or program-agnostic) model "
                "chemistry, or set the method explicitly in the ORCA step."
            )
        if mtype.upper() in ("SQM", "FF", "MLFF"):
            raise RuntimeError(
                f"The model chemistry '{level}' is of type '{mtype}', which ORCA "
                "does not provide. ORCA handles HF, DFT, MP2, and coupled-cluster "
                "(QC) methods."
            )
        if method == "":
            raise RuntimeError(
                f"The model chemistry '{level}' does not name a method ORCA can " "use."
            )
        # If the model chemistry omits a basis, fall back to this node's basis.
        if basis == "":
            basis = P["basis"]
        return method, basis

    def run(self, keywords=None):
        """Run the single-point energy."""
        next_node = super().run(printer)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        printer.important(__(self.description_text(P), indent=self.indent))

        keyword_line = self.keyword_line(P)
        if keywords:
            keyword_line += " " + " ".join(keywords)

        data = self.run_orca(keyword_line)

        self._data = data
        self.analyze(P=P, data=data)

        return next_node

    def analyze(self, indent="", P=None, data=None, **kwargs):
        """Report the energy."""
        if data is None:
            data = getattr(self, "_data", {})
        energy = data.get("energy")
        if energy is None:
            printer.normal(
                __(
                    "ORCA finished but no final single-point energy was found.",
                    indent=self.indent + 4 * " ",
                )
            )
            return

        # Store on the configuration as a property.
        _, configuration = self.get_system_configuration(None)
        try:
            configuration.properties.put(
                "total energy#ORCA#{model}",
                energy,
                _type="float",
                units="E_h",
            )
        except Exception as e:  # pragma: no cover - property store is best-effort
            logger.debug(f"Could not store energy property: {e}")

        printer.normal(
            __(
                f"The total energy is {energy:.8f} E_h.",
                indent=self.indent + 4 * " ",
            )
        )

    def cleanup(self):
        """Nothing to clean up for a single point."""
        pass
