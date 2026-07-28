.. _user-guide:

**********
User Guide
**********

The ORCA plug-in runs `ORCA <https://www.faccts.de/orca/>`_ from a SEAMM
flowchart. It is a *sub-flowchart* plug-in: dropping an **ORCA** step onto the
canvas opens a small flowchart of ORCA sub-steps. Four are available:

* **Energy** — a single-point energy, optionally with the gradient (forces) and
  a range of properties.
* **Optimization** — a geometry optimization (the Energy step plus ORCA's
  ``Opt`` keyword).
* **Frequencies** — the Hessian and harmonic vibrational frequencies, with
  thermochemistry (see below).
* **BSSE** — the counterpoise-corrected energy and gradient of a two-fragment
  complex (see below).

All four share the same controls for the level of theory, so the sections below
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

Initial guess and wavefunction restart
=======================================

Single atoms and open-shell transition metals are often the hardest systems to
converge -- a superposition of atomic densities (ORCA's default ``SAD`` guess)
has little meaning for one atomic center, and can converge to the wrong
electronic state entirely. The **Initial guess** control sets ORCA's SCF
starting guess (``Guess`` in the ``%scf`` block):

* ``default`` -- leave ORCA's own default.
* ``Hueckel``, ``HCore``, ``PAtom``, ``PModel``, ``SAD``, ``SADNO`` -- ORCA's
  built-in guess types. ``PModel`` is often much more stable than ``SAD`` for a
  single atom or ion.
* **Previous wavefunction** -- seed from the orbitals (``orca.gbw``) left by the
  *nearest earlier ORCA step in this flowchart* -- for example, run a cheap,
  robust functional (``PBE``) immediately before a fragile double-hybrid
  (``REVDSD-PBEP86-D4/2021``) on the same atom, and turn this on for the
  double-hybrid step. This only reaches a *different, preceding* node in a
  roughly linear flowchart -- it cannot seed across iterations of a **Loop**,
  since a Loop gives each iteration its own directory (see below).
* **Specified orbitals** -- seed from a named checkpoint file instead, named by
  the **Specified orbitals** control revealed underneath. This is the way to
  seed across loop iterations: for example, a basis-set escalation inside a
  Loop over basis sets, where each iteration reads the checkpoint the
  *previous* iteration (a smaller basis) wrote.

Either wavefunction choice reveals **If wavefunction not found**, controlling
what happens when there is nothing to read (e.g. the very first ORCA step, or
the first pass through a basis-set-escalation loop): ``Throw an error`` (the
default -- an explicit wavefunction request that cannot be honored is treated
as a configuration problem) or one of the ``Use ... guess`` choices, which
falls back to that guess type instead. Set this to a fallback guess to make
"Specified orbitals" safe to leave on for every iteration of a loop.

.. note::

   **ORCA projects across basis sets.** ``MORead`` is not limited to restarting
   with the *same* basis set -- when the checkpoint's basis differs from the
   current job's, ORCA projects the old orbitals onto the new basis (the same
   trick ORCA's own ``Extrapolate`` (CBS) keyword uses internally). So this also
   works for a basis-set escalation, not just same-basis restarts across
   different functionals. It should still be the same atom/system in the same
   electronic state (charge and multiplicity) -- MORead supplies starting
   orbitals, it does not reinterpret the electronic state.

Saving a checkpoint for a later step to read is a separate pair of controls:

* **Save orbital checkpoint** -- ``yes`` copies this run's converged orbitals to
  the file named by **Checkpoint name**, revealed underneath, once the run
  finishes successfully.
* **Checkpoint name** and **Specified orbitals** both resolve the same way:
  ``default`` (the default value of each) derives a label automatically from
  the system's composition, charge, and multiplicity (e.g. ``Co_q0_m4``) -- the
  usual choice, since it naturally gives one checkpoint per atom/electronic
  state, shared across an inner loop (e.g. over basis sets) for that atom, and
  reset automatically when an outer loop moves to a new atom. Any other bare
  name is stored under a ``checkpoints`` folder at the top of the job; an
  absolute path is used as-is, e.g. to keep a checkpoint outside this job and
  reuse it across separate flowchart runs.

.. note::

   **Reading a checkpoint from another job.** Only **Specified orbitals** can
   reference *another* job -- ``job://<job number>/<name>`` -- since a job
   must never write into another job (so this is not available for
   **Checkpoint name**). Use ``job://<job number>/default`` when you want that
   other job's own auto-derived name for this same system but do not know
   what it resolved to -- for example, ``job://53/default`` picks up whatever
   ``Co_q0_m4``-style label job 53 auto-derived for this atom. This needs
   SEAMM's managed ``Jobs/<project>/Job_NNNNNN`` directory layout (i.e. jobs
   run through the dashboard/JobServer, not a bare local run) to locate the
   other job.

For a basis-set-escalation loop (e.g. looping ``def2-SVP`` -> ``def2-TZVP`` ->
``def2-QZVP`` for one atom), turn on both **Save orbital checkpoint** and
**Initial guess** = **Specified orbitals** (leaving **Checkpoint name** /
**Specified orbitals** on ``default`` so they agree automatically), and set
**If wavefunction not found** to a fallback guess so the first iteration -- which
has nothing to read yet -- does not error. Loop from the smallest basis to the
largest: projecting a converged small-basis density up into a larger virtual
space is cheap and well-conditioned, which is the direction ORCA's own CBS
extrapolation uses internally.

Extra keywords and extra ORCA blocks
=====================================

Two free-text escape hatches cover anything without a dedicated control:

* **Extra keywords** -- additional ORCA ``!`` keywords, appended after
  everything the other controls generate (e.g. ``RIJCOSX``, ``NoFrozenCore``,
  ``SlowConv``).
* **Extra ORCA blocks** -- literal ORCA input (one or more ``%`` blocks),
  inserted verbatim right before the geometry, after the blocks the controls
  above generate. Use this for one-off SCF-stabilization tricks with no
  dedicated control, for example on a difficult atom::

     %scf
       MaxIter 400
       Shift Shift 0.3 ErrStart 0.05 end
       DIISBfac 1.1
     end

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

Geometry optimization
=====================

The **Optimization** sub-step adds ORCA's ``Opt`` keyword to the energy run and,
when it finishes, stores the optimized geometry. A **Convergence** control
selects ORCA's optimization preset (``LooseOpt`` … ``VeryTightOpt``).

A **Structure handling** control chooses where the optimized geometry goes. It
defaults to *Overwrite the current configuration*, so the optimized structure
flows straight into any following sub-step — for example a **Frequencies**
calculation runs at the optimized geometry without any extra wiring. You can
instead create a new configuration (or a new system and configuration) to keep
the starting structure alongside the optimized one, or discard the optimized
structure entirely.

Frequencies (Hessian and thermochemistry)
==========================================

The **Frequencies** sub-step computes the Hessian and the harmonic vibrational
frequencies, the IR intensities, and the thermochemistry (zero-point energy,
thermal enthalpy, entropy, and Gibbs free energy). The level-of-theory controls
are the same as the Energy step, plus two extra controls:

* **Second derivatives** — ``analytic`` (ORCA's ``AnFreq``; much faster, needs
  an analytic second derivative for the method — HF, most DFT functionals, MP2)
  or ``numerical`` (``NumFreq``; finite-differences the gradient, works for any
  method with a gradient but is considerably more expensive).
* **Temperature** — the temperature for the thermochemistry (the pressure is
  ORCA's default of 1 atm).
* **Structure handling** — where to store the structure and its properties. A
  frequency calculation does not change the geometry, so this defaults to
  *Overwrite the current configuration*; you can instead create a new
  configuration (or a new system and configuration) to hold the results, or
  discard the structure.

Tick the results you want on the Results tab: the **frequencies**, the **IR
intensities**, the **zero-point energy**, the **enthalpy**, the **Gibbs free
energy**, the **number of imaginary frequencies**, and the **largest zero-mode
frequency**. Imaginary (negative) frequencies are reported and flagged — a
minimum has none, a transition state has one.

The output lists the vibrational frequencies and their IR intensities in a
table. The frequencies (and IR intensities) are also written to
``frequencies.csv`` in the step's directory for easy access, and an
``IR_spectrum.graph`` file is written with the IR spectrum as a stick trace plus
a Lorentzian-broadened trace (FWHM ≈ 15 cm⁻¹) that mimics an experimental
spectrum. Open the ``.graph`` file in the SEAMM dashboard to view it
interactively.

.. note::

   **Units.** Following SEAMM's SI-based convention, the zero-point energy,
   enthalpy, and Gibbs free energy are reported in **kJ/mol**. (Orbital energies
   — HOMO/LUMO — are reported in **eV** throughout ORCA, the conventional unit
   for orbital energies.)

The report also lists the **largest of the 5 or 6 nominally-zero
translation/rotation frequencies**. These should be zero; how far the largest
one departs from zero is a convenient gauge of the numerical accuracy of the
Hessian. ORCA *projects* the translations and rotations out of its printed
frequencies (it reports them as exactly ``0.00``), so this residual is instead
obtained by diagonalizing the **raw, un-projected** mass-weighted Cartesian
Hessian from ``orca.hess`` — the value you see is therefore the true numerical
residual, not zero.

.. note::

   The same analytic second derivative also backs the ORCA **MDI engine**: it
   answers a custom ``<HESSIAN`` command (``AnFreq``) so a driver — e.g. the
   Normal Mode Sampling step — can pull the analytic Hessian over a warm MDI
   connection. The engine advertises ``<HESSIAN`` **only when ORCA has an
   analytic Hessian for the method** (HF, MP2, ordinary DFT) and not for double
   hybrids or ``(DLPNO-)CCSD(T)``, so the driver's capability check is truthful:
   it uses the analytic Hessian when offered and finite-differences the forces
   otherwise.

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
