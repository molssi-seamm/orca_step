.. _user-guide:

**********
User Guide
**********

The ORCA plug-in runs `ORCA <https://www.faccts.de/orca/>`_ from a SEAMM
flowchart. It is a *sub-flowchart* plug-in: dropping an **ORCA** step onto the
canvas opens a small flowchart of ORCA sub-steps. Three are available:

* **Energy** — a single-point energy, optionally with the gradient (forces) and
  a range of properties.
* **Optimization** — a geometry optimization (the Energy step plus ORCA's
  ``Opt`` keyword).
* **BSSE** — the counterpoise-corrected energy and gradient of a two-fragment
  complex (see below).

All three share the same controls for the level of theory, so the sections below
apply to any of them.

Choosing the level of theory
============================

Each sub-step can take its method and basis set either from a preceding
**Model Chemistry** step (the default — leave *Use the global model chemistry*
on) or set them explicitly. Turn *Use the global model chemistry* off to choose
them in the step itself.

Method
------

The **Method** pull-down offers Hartree–Fock, MP2/RI-MP2, coupled cluster
(``CCSD(T)``, ``DLPNO-CCSD(T)``), the explicitly-correlated **F12** variants
(``CCSD(T)-F12D/RI``, ``DLPNO-CCSD(T)-F12D``), and **DFT**. For everything except
DFT the method *is* the ORCA keyword. Choosing **DFT** reveals two further,
indented controls:

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

Explicitly-correlated (F12) methods
-----------------------------------

The F12 methods (``CCSD(T)-F12D``, ``DLPNO-CCSD(T)-F12D``) reach near
complete-basis-set accuracy from a modest basis, but they **must** be paired
with one of ORCA's F12-optimized orbital bases — ``cc-pVDZ-F12``,
``cc-pVTZ-F12``, or ``cc-pVQZ-F12`` — and each of those needs a matching
*complementary auxiliary basis set* (CABS). The step **adds the CABS
automatically**: it derives ``<basis>-CABS`` from the F12 basis you choose (so
``cc-pVTZ-F12`` → ``cc-pVTZ-F12-CABS``, and likewise for DZ/QZ), unless you have
already put a CABS in the extra keywords. Leave the auxiliary basis on
``AutoAux`` — it supplies the remaining RI fitting bases; the CABS is separate
and is not something AutoAux can generate.

To keep this foolproof, **selecting an F12 method narrows the basis-set list to
just the F12 bases** (and forces the ORCA-internal source), so you can only pick
a valid one; switching back to a non-F12 method restores the full list. If a
hand-edited flowchart still pairs an F12 method with a non-F12 basis (and no
CABS in the extra keywords), the run is stopped early with a clear message.

.. tip::

   F12 methods largely **eliminate basis-set superposition error**, so the
   counterpoise (BSSE) correction is essentially zero for them (often a few
   thousandths of a kcal/mol — at the numerical-noise level). With F12 you can
   skip the **BSSE** sub-step and take the interaction energy directly from a
   plain **Energy** run, which is ~5× cheaper. Reserve the BSSE step for
   conventional finite-basis methods (HF, DFT, MP2, canonical CCSD(T)), where the
   correction is real.

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

Integration grid
================

The **Integration grid** control sets ORCA's numerical-integration grid preset
(used for the DFT exchange-correlation integrals and the RIJCOSX/COSX grids):

* ``default`` — leave ORCA's own default (``DEFGRID2``), which is robust for most
  work;
* ``DEFGRID1`` — coarser and faster (roughly the old ORCA-4 default accuracy);
* ``DEFGRID2`` — the current default;
* ``DEFGRID3`` — finer and more conservative, for cases sensitive to the grid
  (e.g. dispersion-dominated energies, tight convergence, or numerical noise in
  forces).

It only affects methods that use a grid (DFT, RIJCOSX) and is ignored otherwise.
For a machine-learned-force-field training set, a finer grid (``DEFGRID3``) can
be worth the cost to keep the energy/force surface smooth.

When the control is left on ``default``, the step **automatically switches to
``DEFGRID3`` for basis sets with high angular momentum** — h functions or above,
such as ``cc-pV5Z`` — which the default ``DEFGRID2`` integrates less accurately.
The angular momentum is read from the Basis Set Exchange (which also covers
ORCA's internal basis names), so it works for both internal and BSE bases;
choosing a grid explicitly always overrides this.

SCF convergence
===============

The **SCF convergence** control sets ORCA's SCF convergence-tolerance preset,
from loosest to tightest: ``SLOPPYSCF`` → ``LOOSESCF`` → ``NORMALSCF`` →
``STRONGSCF`` → ``TIGHTSCF`` → ``VERYTIGHTSCF`` → ``EXTREMESCF``. ``default``
leaves ORCA's own default (``NORMALSCF`` for a single point; ORCA already
tightens it to ``TIGHTSCF`` for optimizations).

The step defaults to **``TIGHTSCF``**, which is a good choice for smooth
energies and forces (and matches what ORCA uses for optimizations). Loosen it
only to save time when high precision is not needed; tighten it (``VERYTIGHTSCF``
/ ``EXTREMESCF``) for very-high-accuracy or numerically delicate work.

SCF SThresh (linear dependence)
===============================

The **SCF SThresh** control sets ORCA's ``%scf SThresh`` value — the threshold
below which an eigenvalue of the overlap matrix is treated as zero, so the
corresponding (near-)linearly-dependent combination of basis functions is
dropped from the calculation. Redundant, near-linearly-dependent basis functions
cause numerical instabilities in the SCF, and diffuse-heavy sets (the augmented
``aug-cc-pVXZ`` / ``ma-…`` families in particular) are the usual culprits.
Raising ``SThresh`` above the default removes more of these functions and can
cure the resulting SCF trouble.

* ``default`` — emit nothing, so the SCF-convergence preset above (or ORCA's own
  default of ``1.0e-07``) governs ``SThresh``.
* an explicit value — written to ``%scf SThresh``, overriding whatever the preset
  would otherwise set. ORCA recommends keeping it between ``1e-8`` and ``1e-5``;
  raise it toward ``1e-6`` to shed near-dependent functions when a diffuse basis
  will not converge.

.. caution::

   Values beyond ``1e-6`` must be used carefully in **geometry optimizations**
   and when **comparing conformers**: because different basis functions can be
   cut off at different geometries, the final basis set — and hence the energy —
   can vary discontinuously from one structure to the next.

See the ORCA manual's `basis-set section
<https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/basisset.html>`_
for the full discussion of linear dependence and its automatic removal.

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

Counterpoise (BSSE) corrections
===============================

The **BSSE** sub-step computes the counterpoise-corrected (Boys--Bernardi)
energy **and gradient** of a complex of two fragments — removing the basis-set
superposition error that artificially over-stabilizes the interaction — in a
single ORCA run. It is aimed at machine-learned-force-field (MLFF) training
data that is BSSE-free on both the energy surface and the forces. Internally it
drives ORCA's *Compound* facility (the ``BSSEGradient`` script by D. G. Liakos &
F. Neese), which runs the five sub-calculations (the dimer, and each fragment
both in the full dimer basis and in its own basis) and assembles the correction.

Use it with the conventional finite-basis methods (HF, DFT, MP2, canonical
CCSD(T)), where BSSE is a real effect. For **explicitly-correlated F12 methods**
the superposition error is already negligible, so the counterpoise correction is
essentially zero and this step is unnecessary — take the interaction energy
directly from a plain **Energy** run instead (see *Explicitly-correlated (F12)
methods* above).

The level-of-theory controls are the same as the Energy step. Three extra
controls define the correction:

* **Fragments** — how to split the complex into the two fragments. ``auto (2
  molecules)`` (the default) uses the two separate molecules in the structure
  (so it works directly on a dimer from the **Dimer Builder** step), and errors
  if there are not exactly two. ``specified`` takes fragment A from the
  **Fragment A atoms** field, given as **atom numbers** (1-based, as shown in the
  structure) — a comma/space separated list and/or ranges, e.g. ``1-3, 5``; the
  remaining atoms form fragment B. The **Fragment A atoms** field is shown only
  when ``specified`` is selected.
* **Compute the gradient** — ``yes`` (the default) computes the
  counterpoise-corrected gradient (forces) as well as the energy, which is what
  MLFF training needs. ``no`` (energy only) is cheaper and, importantly, allows
  methods that have **no analytic gradient** in ORCA — notably ``CCSD(T)`` and
  ``DLPNO-CCSD(T)`` — for gold-standard counterpoise *interaction energies*.
* **Optimize free monomers** — whether to relax each isolated monomer before
  taking the correction. Leave it ``no`` for a fixed-geometry PES / MLFF target.
* **Write the wavefunction (wfx) file** — default ``no``. When ``yes``, the
  dimer's density is kept and converted (via ``orca_2aim``) to an ``orca.wfx``,
  so a following **Atomic Charges** step can partition it into DDEC6 charges on
  the CP complex — exactly as after an Energy step.

The step reports, and offers on the Results tab, the **BSSE-corrected energy**
(the ``energy`` result), the **uncorrected** (raw) complex energy, the **BSSE
correction** (corrected minus uncorrected, also shown in kcal/mol), and — when
requested — the corrected **gradient**. Tick any of them on the Results tab to
save it to a variable, table, or the property database.

.. note::

   **Phase-1 limitations.** This first version supports a **neutral,
   closed-shell** complex of **exactly two** fragments, using an
   **ORCA-internal** basis set (not the Basis Set Exchange). For the **energy**,
   any method works — HF, DFT (including dispersion-corrected and
   **double-hybrid** functionals such as ``REVDSD-PBEP86-D4/2021``), MP2, and
   ``(DLPNO-)CCSD(T)``. Computing the **gradient** additionally requires an
   **analytic** gradient, so the numerical-gradient methods (``(DLPNO-)CCSD(T)``)
   are available in *energy-only* mode only. The same charge/multiplicity is
   applied to each monomer, so charged or open-shell fragments are refused with a
   clear message; the SThresh control and the extra property analyses (bond
   orders, Hirshfeld charges, polarizability, saved wavefunction) are not
   available in the BSSE sub-step. N-fragment and code-agnostic counterpoise
   corrections are planned as a separate, general BSSE step.

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

By default ORCA uses all the cores the machine or batch job provides. Settings
come from two files with different jobs:

* **How to run ORCA** lives in ``~/SEAMM/orca.ini`` -- the full path to the
  executable and the OpenMPI ``library-path`` for parallel runs:

  .. code-block:: ini

     [local]
     installation = local
     code = /path/to/orca
     library-path = /path/to/orca-mpi/lib

* **User run options** live in the ``[orca-step]`` section of the main SEAMM
  configuration (``~/.seamm.d/seamm.ini``), and can also be given on the command
  line:

  .. code-block:: ini

     [orca-step]
     ncores = available        # or an integer, or 1 to force serial
     memory = available        # or 'all', or e.g. '3 GB' (per process)

* **ncores** — how many processes ORCA may use (its ``%pal``). ``available``
  (the default) uses all cores the machine/job provides; give an integer to cap
  it, or ``1`` to force serial.
* **memory** — the per-process memory for ORCA's ``%maxcore``. ``available``
  (the default) scales to the memory per core; ``all`` divides the whole node
  among the processes; or give an explicit amount such as ``3 GB``.
* **library-path** (in ``orca.ini``) — the ``lib`` directory of the OpenMPI that
  **matches the version ORCA was built against**. ORCA 6.1 requires OpenMPI
  4.1.x and does *not* support 5.x; mixing versions makes parallel runs abort
  with a
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

Driving ORCA as an MDI engine
=============================

Besides running as an ordinary flowchart step, ORCA can act as a persistent
`MDI <https://molssi-mdi.github.io/MDI_Library/>`_ engine for steps that set up a
*model chemistry* and then evaluate it at many geometries -- for example the
**Dimer Builder** step's energy-based contact search. You do not configure this
in the ORCA step itself: put a **Model Chemistry** step in the flowchart, choose
an ORCA model chemistry there (e.g. ``ORCA:DFT@B3LYP/def2-SVP``), and the driving
step launches ORCA as the engine automatically.

Because ORCA has no in-process interface, the engine runs the ``orca`` binary
once per geometry in a persistent working directory, **reusing the previous
geometry's orbitals** (``orca.gbw``) as the SCF guess -- the main saving for a
series of nearby structures. Two consequences:

* **Only methods with an analytic gradient are offered via MDI.** The engine
  always computes the energy and forces together (``EnGrad``), so
  ``DLPNO-CCSD(T)``, ``CCSD(T)`` and the non-self-consistent ``wB97M(2)`` /
  ``wB97X-2`` are *not* MDI-capable; choosing one gives a clear error. Everything
  else (HF, MP2, and the analytic-gradient functionals) works.
* **ORCA is not start-up-dominated**, so the per-geometry cost is real -- the MDI
  benefit here is orbital reuse and a uniform interface, not the large speed-up
  that cheap engines see. Use an inexpensive functional for jobs that only need
  the energy surface to guide them (such as contact finding); reserve expensive,
  high-accuracy calculations for ordinary single-point steps.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
