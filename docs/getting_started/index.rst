***************
Getting Started
***************

Installation
============
The ORCA step is installed with the `SEAMM Manager`_, and is probably already part of
your SEAMM installation. To add it, or bring it up to date::

  seamm-manager install orca-step
  seamm-manager update orca-step

or use the Manager's window. ORCA itself is licensed software that you install
yourself, so installing the step does not install ORCA. Instead it creates
``~/SEAMM/orca.ini`` -- or ``orca.ini`` in whichever SEAMM installation you are
working on -- which tells SEAMM where ORCA is and how to run it. An existing file
is never changed.

.. _SEAMM Manager: https://molssi-seamm.github.io/getting_started/installation/seamm-manager.html

If ``code`` is not given in ``orca.ini``, SEAMM runs the ``orca`` on your ``PATH``,
using its full path as ORCA requires.

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
