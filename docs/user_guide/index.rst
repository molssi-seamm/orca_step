.. _user-guide:

**********
User Guide
**********

The ORCA plug-in runs `ORCA <https://www.faccts.de/orca/>`_ from a SEAMM
flowchart. It is a *sub-flowchart* plug-in: dropping an **ORCA** step onto the
canvas opens a small flowchart of ORCA sub-steps. Two are available:

* **Energy** — a single-point energy, optionally with the gradient (forces) and
  a range of properties.
* **Optimization** — a geometry optimization (the Energy step plus ORCA's
  ``Opt`` keyword).

Both share the same controls for the level of theory, so the sections below
apply to either.

Choosing the level of theory
============================

Each sub-step can take its method and basis set either from a preceding
**Model Chemistry** step (the default — leave *Use the global model chemistry*
on) or set them explicitly. Turn *Use the global model chemistry* off to choose
them in the step itself.

Method
------

The **Method** pull-down offers Hartree–Fock, MP2/RI-MP2, coupled cluster
(``CCSD(T)``, ``DLPNO-CCSD(T)``), and **DFT**. For everything except DFT the
method *is* the ORCA keyword. Choosing **DFT** reveals two further, indented
controls:

* **Functional type** — ORCA's own classification of the functional: local,
  GGA, meta-GGA, (range-separated) hybrid, and (range-separated) double-hybrid.
* **Functional** — the functional itself, filtered to the chosen type.

Every functional ORCA documents is available, grouped by type; picking a type
narrows the functional list so it stays readable. Double-hybrid functionals
(for example ``REVDSD-PBEP86-D4/2021`` or ``B2PLYP``) need an auxiliary ``/C``
fitting basis for their MP2 part — leave the *Auxiliary (fitting) basis* on
``AutoAux`` and it is supplied automatically.

Basis set
---------

Type any basis name ORCA knows internally, pick one of the curated families
from the list (Pople, Dunning correlation-consistent, and Karlsruhe def2,
including the diffuse ``…D`` and minimally-augmented ``ma-…`` variants), or
press **...** to browse the `Basis Set Exchange
<https://www.basissetexchange.org/>`_ on a periodic table, filtered to the
elements in your system. Selecting **Basis Set Exchange** as the *Basis set
source* opens the same picker. A basis chosen from the Exchange is stored as
``bse:NAME`` and embedded in the ORCA input, so the definition is identical
across codes. The list is ordered by family and, within a family, into
valence / polarization / diffuse ladders that each rise DZ → TZ → QZ → 5Z, so a
sensible progression is a single ladder read top to bottom.

Complete-basis-set (CBS) extrapolation
--------------------------------------

Set **Basis-set extrapolation (CBS)** to ``2/3``, ``3/4``, or ``4/5`` to
extrapolate the energy to the complete-basis-set limit from two successive
cardinal numbers (double/triple, triple/quadruple, quadruple/quintuple), using
the **Extrapolation family** (``cc``, ``aug-cc``, ``def2``, or ``ANO``). This
is a *single* ORCA job (its ``Extrapolate`` keyword): ORCA runs both basis sets
and extrapolates the SCF and correlation contributions separately. When it is
on, the fixed basis set above is ignored. Note that an extrapolated energy has
**no gradient** in ORCA, so extrapolation and the *gradients* result are
mutually exclusive.

Energies, gradients, and forces
===============================

The **Results** tab lists everything the step can produce. Tick a result to
save it to a variable, a table, or the structure's property database (see
below).

Requesting **gradients** makes the step compute the nuclear gradient (i.e. the
forces). ORCA computes an *analytic* gradient (``EnGrad``) when one exists for
the chosen method, and automatically falls back to a *numerical* gradient
(``NumGrad``) when it does not — for example ``DLPNO-CCSD(T)`` (no analytic
``(T)`` gradient) and the non-self-consistent ``wB97M(2)`` / ``wB97X-2`` double
hybrids. A note is printed when the (more expensive) numerical gradient is used.
Most functionals — including the standard double hybrids — have analytic
gradients, so a single-point DFT run yields the energy **and** forces cheaply,
which is convenient for generating machine-learned force-field training data.

Saving results to the database
==============================

Both scalar results (the energies, HOMO/LUMO energies and gap, the dipole
magnitude, ``<S^2>``, the polarizability) and array results (the **gradient**,
the dipole-moment vector, the Mulliken/Löwdin/Hirshfeld charges, the Mayer
valences, and the rotational constants) can be stored as properties on the
configuration. In the Results tab, tick the database column for the result; it
is stored under a name like ``gradients#ORCA#<model>``, where ``<model>`` is the
level of theory.

Other properties and outputs
============================

The Energy step can also compute Mayer bond orders and Hirshfeld charges (and
optionally apply them to the structure), the dipole polarizability, and an
analytic wavefunction (``.wfx``, via ``orca_2aim``) for a following **Atomic
Charges** step to partition into DDEC6 charges. Every run cites ORCA, the DFT
functional, the basis set, and the supporting integral / exchange-correlation
libraries.

Running ORCA in parallel
========================

By default ORCA now uses all the cores the machine or batch job provides. The
step reads its resource settings from the ``[orca-step]`` section of
``~/SEAMM/orca.ini`` (and from command-line options of the same name):

.. code-block:: ini

   [local]
   code = /path/to/orca

   [orca-step]
   ncores = available        # or an integer, or 1 to force serial
   memory = available        # or 'all', or e.g. '3 GB' (per process)
   library-path = /path/to/orca/openmpi/lib

* **ncores** — how many processes ORCA may use (its ``%pal``). ``available``
  (the default) uses all cores the machine/job provides; give an integer to cap
  it, or ``1`` to force serial.
* **memory** — the per-process memory for ORCA's ``%maxcore``. ``available``
  (the default) scales to the memory per core; ``all`` divides the whole node
  among the processes; or give an explicit amount such as ``3 GB``.
* **library-path** — the ``lib`` directory of the OpenMPI that **matches the
  version ORCA was built against**. ORCA 6.1 requires OpenMPI 4.1.x and does
  *not* support 5.x; mixing versions makes parallel runs abort with a
  ``BLAS-ERROR``. A dedicated conda env is the easy way to get the right one::

     conda create -n orca-mpi -c conda-forge "openmpi=4.1"

  then set ``library-path`` to that env's ``lib`` directory. The step
  automatically puts the matching ``mpirun`` (the sibling ``bin`` directory) on
  ``PATH`` so ORCA launches its workers with the correct OpenMPI.

.. note::

   **macOS.** ORCA does not pass the ``DYLD_*`` loader variables to the MPI
   processes it spawns, so ``library-path`` alone does not let the dynamic
   loader find ``libmpi`` on a Mac. Put the OpenMPI libraries on the default
   search path once, for example by symlinking into ``/usr/local/lib``::

      ln -s /path/to/orca-mpi/lib/libmpi.40.dylib /usr/local/lib/

   On Linux the exported ``LD_LIBRARY_PATH`` is inherited normally, so
   ``library-path`` is sufficient there.

If you do not have a matching OpenMPI, set ``ncores = 1`` to run serially.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
