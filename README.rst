==================
SEAMM ORCA Plug-in
==================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/orca_step
   :target: https://github.com/molssi-seamm/orca_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/orca_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/orca_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/orca_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/orca_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/orca_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/orca_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/orca_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/orca_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/orca_step.svg
   :target: https://pypi.python.org/pypi/orca_step
   :alt: PyPi VERSION

A SEAMM plug-in for ORCA

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/orca_step/index.html
* Code: https://github.com/molssi-seamm/orca_step

Features
--------

A SEAMM plug-in for ORCA, a general-purpose quantum-chemistry program, with an
emphasis on accurate molecular calculations such as DLPNO-CCSD(T).

Like the MOPAC and Gaussian steps, the ORCA step is a sub-flowchart: you add an
ORCA node to your flowchart and then build a small sub-flowchart of ORCA
capabilities inside it. Initially the available capabilities are:

* **Energy** -- a single-point energy.
* **Optimization** -- a geometry optimization.

Further capabilities (frequencies, properties, ...) will be added.

Methods are described by the step's metadata and can be set either explicitly in
the ORCA dialog (similar to the Gaussian step) or, by default, taken from a
preceding **Model Chemistry** step. Basis sets default to ORCA's built-in
families (Pople, Dunning ``cc``, and Karlsruhe ``def2``), with the Basis Set
Exchange available as a planned opt-in source.

Acknowledgements
----------------

This package was created with the `molssi-seamm/cookiecutter-seamm-plugin`_ tool, which
is based on the excellent Cookiecutter_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award CHE-2136142.

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
