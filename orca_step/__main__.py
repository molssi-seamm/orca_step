# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""Handle the installation of the ORCA step."""

from .installer import Installer


def run():
    """Find and/or register the ORCA executable in seamm.ini."""
    installer = Installer()
    installer.run()


if __name__ == "__main__":
    run()
