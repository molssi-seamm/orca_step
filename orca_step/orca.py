# -*- coding: utf-8 -*-

"""The main ORCA node: holds and drives a sub-flowchart of ORCA capabilities."""

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


class ORCA(orca_step.ORCABase):
    """The main ORCA node. Like the MOPAC and Gaussian nodes, it owns a
    sub-flowchart whose nodes are ORCA capabilities (Energy, Optimization, ...).

    See Also
    --------
    TkORCA, ORCABase
    """

    def __init__(
        self,
        flowchart=None,
        title="ORCA",
        namespace="org.molssi.seamm.orca",
        extension=None,
        logger=logger,
    ):
        """Initialize the ORCA node and its sub-flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The flowchart that contains this step.
        title: str
            The name displayed in the flowchart.
        namespace: str
            The stevedore namespace for the sub-flowchart's plug-ins.
        extension: None
            Not yet implemented.
        """
        logger.debug(f"Creating ORCA {self}")

        self.subflowchart = seamm.Flowchart(
            parent=self, name="ORCA", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = orca_step.metadata

    def set_id(self, node_id):
        """Set the id of this node and propagate to the sub-flowchart."""
        next_node = super().set_id(node_id)
        self.subflowchart.set_ids(self._id)
        return next_node

    def description_text(self, P=None):
        """Describe what the ORCA sub-flowchart will do."""
        self.subflowchart.root_directory = self.flowchart.root_directory

        node = self.subflowchart.get_node("1").next()
        text = self.header + "\n\n"
        while node is not None:
            try:
                text += __(node.description_text(), indent=3 * " ").__str__()
            except Exception as e:
                logger.warning(f"Error describing ORCA flowchart: {e} in {node}")
            text += "\n"
            node = node.next()
        return text

    def run(self):
        """Run the ORCA sub-flowchart, node by node."""
        # The node after this one, to return at the end
        next_node = super().run(printer)

        printer.important(self.header)
        printer.important("")

        self.subflowchart.root_directory = self.flowchart.root_directory
        self.subflowchart.executor = self.flowchart.executor

        node = self.subflowchart.get_node("1").next()
        while node is not None:
            if node.is_runable:
                node.run()
                # Separate each sub-step's output with a blank line.
                printer.normal("")
            node = node.next()

        return next_node
