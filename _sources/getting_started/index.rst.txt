***************
Getting Started
***************

Installation
============
The ORCA step is probably already installed in your SEAMM
environment, but if not or if you wish to check, follow the directions for the `SEAMM
Installer`_. The graphical installer is the easiest to use. In the SEAMM conda
environment, simply type:: 

  seamm-installer

or use the shortcut if you installed one. Switch to the second tab, `Components`, and
check for `orca-step`. If it is not installed, or
can be updated, check the box next to it and click `Install selected` or `Update
selected` as appropriate.

The non-graphical installer is also straightforward::

  seamm-installer install --update orca-step

will ensure both that it is installed and up-to-date.

.. _SEAMM Installer: https://molssi-seamm.github.io/installation/index.html

A first calculation
===================
Add an **ORCA** step to your flowchart and open it to reveal the ORCA
sub-flowchart, then add an **Energy** (or **Optimization**) sub-step. In that
sub-step either leave *Use the global model chemistry* on (to take the method
and basis from a preceding **Model Chemistry** step) or turn it off and choose
them directly — for example set **Method** to ``DFT``, **Functional type** to
*global hybrid*, **Functional** to ``B3LYP``, and the basis to ``def2-TZVP``.
On the **Results** tab, tick the energy (and *gradients* if you want forces) to
save them. Run the flowchart as usual.

To run ORCA in parallel, set ``ncores`` in the ``[orca-step]`` section of the
main SEAMM configuration (``~/.seamm.d/seamm.ini``); the path to ORCA and, for
parallel runs, the OpenMPI ``library-path`` go in ``~/SEAMM/orca.ini`` (see the
User Guide).

That should be enough to get started. For more detail about the functionality in this plug-in, see the :ref:`User Guide <user-guide>`.
