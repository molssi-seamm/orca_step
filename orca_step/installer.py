# -*- coding: utf-8 -*-

"""Installer for the ORCA plug-in.

ORCA is licensed software installed manually from https://www.faccts.de / the
ORCA forum, so this installer locates the existing ``orca`` executable and
registers its full path in seamm.ini (ORCA needs its full path to launch its
sub-programs).
"""

import importlib
import logging
import shutil

import seamm_installer

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Locate the ORCA executable and register it in seamm.ini."""

    def __init__(self, logger=logger):
        super().__init__(logger=logger)
        logger.debug("Initializing the ORCA installer object.")

        self.section = "orca-step"
        self.executables = ["orca"]
        self.resource_path = importlib.resources.files("orca_step") / "data"

    def exe_version(self, config):
        """Return the ORCA name and version.

        ORCA has no stable ``--version`` flag (the version is in the run banner),
        so this just confirms the executable is reachable.
        """
        code = config.get("code", "orca")
        path = shutil.which(code) or code
        return "ORCA", "unknown" if path is None else "unknown"
