# -*- coding: utf-8 -*-

"""
orca_step
A SEAMM plug-in for ORCA, with a sub-flowchart of capabilities (Energy,
Optimization, ...), modeled on the MOPAC and Gaussian steps.
"""

# Metadata and the shared base first (others import them).
from .metadata import metadata  # noqa: F401
from .orca_base import ORCABase  # noqa: F401

# Sub-step parameters / nodes / Tk / factories.
from .energy_parameters import EnergyParameters  # noqa: F401
from .energy import Energy  # noqa: F401
from .tk_energy import TkEnergy  # noqa: F401
from .energy_step import EnergyStep  # noqa: F401

from .optimization_parameters import OptimizationParameters  # noqa: F401
from .optimization import Optimization  # noqa: F401
from .tk_optimization import TkOptimization  # noqa: F401
from .optimization_step import OptimizationStep  # noqa: F401

from .bsse_parameters import BSSEParameters  # noqa: F401
from .bsse import BSSE  # noqa: F401
from .tk_bsse import TkBSSE  # noqa: F401
from .bsse_step import BSSEStep  # noqa: F401

# Main node, Tk, and factory.
from .orca import ORCA  # noqa: F401
from .tk_orca import TkORCA  # noqa: F401
from .orca_step import ORCAStep  # noqa: F401
from .orca_step import mc_method_alias, mc_method_unalias  # noqa: F401

# Versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
