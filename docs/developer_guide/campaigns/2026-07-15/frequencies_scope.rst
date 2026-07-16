Frequencies sub-step and <HESSIAN -- scope and design
=====================================================

Goal
----

Two related additions to the ORCA plug-in:

1. A **Frequencies** sub-step (a subclass of Energy) that runs ORCA's ``Freq``
   (analytic ``AnFreq``) or ``NumFreq`` (numerical) second derivatives and
   reports the Hessian, harmonic vibrational frequencies, IR intensities, and
   the thermochemistry (zero-point energy, enthalpy, entropy, Gibbs free energy)
   at a chosen temperature. Imaginary frequencies are reported and flagged.

2. A custom **``<HESSIAN``** MDI command on the ORCA MDI engine that returns the
   analytic Cartesian Hessian, so a driver holding a *warm* engine can pull the
   Hessian without a fresh process launch.

The immediate driver is the ``normal_mode_sampling_step`` plug-in, which samples
the normal-mode space of a molecule and needs the Hessian at the optimized
monomer. It obtains that Hessian through ``seamm_mdi.MDIEngine``: an analytic
Hessian over ``<HESSIAN`` when the engine offers one, otherwise a
finite-difference of the forces over the same warm engine.

The Frequencies sub-step
------------------------

``Freq``/``NumFreq`` appended to the Energy input; ORCA writes the mass-weighted
Hessian to ``orca.hess`` and the vibrational analysis to ``orca.out``. Beyond the
frequencies and IR intensities, the step reports the **largest of the 5 or 6
nominally-zero translation/rotation frequencies** as a gauge of the numerical
quality of the Hessian -- computed from the raw, un-projected ``orca.hess`` (ORCA
prints those modes as exactly ``0.00``, which hides the real residual). It writes
``frequencies.csv`` and an ``IR_spectrum.graph`` (stick + Lorentzian-broadened
trace), and honours the standard Structure-handling options.

The <HESSIAN MDI command
------------------------

The bundled engine ``data/orca_mdi.py`` answers ``<HESSIAN`` by running an
``! AnFreq`` job for the current geometry, parsing ``orca.hess``, and sending the
**3N x 3N Cartesian Hessian in hartree/bohr^2, row-major**. The result is cached
and invalidated when the geometry, charge, or multiplicity changes, so repeated
requests at one point are free.

Truthful capability advertising (the key design decision)
---------------------------------------------------------

MDI has no standard Hessian node, so a driver must *discover* the capability at
run time (``MDI_Check_Command_Exists`` / ``MDIEngine.supports("<HESSIAN")``). For
that discovery to be trustworthy, **the engine advertises ``<HESSIAN`` only when
ORCA genuinely has an analytic Hessian for the method**:

* ``get_mdi_engine_command`` passes ``--hessian yes/no`` from
  :func:`method_has_analytic_hessian`;
* HF, MP2, and the ordinary DFT functionals have an analytic Hessian -> advertise;
* **double hybrids** (e.g. revDSD-PBEP86-D4) have an analytic *gradient* but **no
  analytic Hessian** (they need ``NumFreq``) -> do **not** advertise;
* **(DLPNO-)CCSD(T)** and the non-self-consistent double hybrids have neither ->
  do not advertise.

An earlier iteration advertised ``<HESSIAN`` unconditionally and let ``<HESSIAN``
degrade to a numerical Hessian when no analytic one existed. That was a bad
design: the capability check lied, and a driver expecting a cheap analytic
Hessian could unknowingly trigger a ruinously expensive ``NumFreq`` (catastrophic
for coupled cluster). The engine now **never silently falls back to NumFreq over
MDI** -- if it does not advertise ``<HESSIAN``, the driver sees ``supports()``
return ``False`` and finite-differences the forces itself (its own explicit,
visible choice). Use the Frequencies *sub-step* when a numerical Hessian is
actually wanted.

Cross-repo pieces (all released 2026-07-15/16)
----------------------------------------------

* **orca_step 2026.7.15** -- Frequencies sub-step + the ``<HESSIAN`` engine
  command (this campaign); **2026.7.16** -- structure handling, kJ/mol units,
  the zero-mode residual, the IR-spectrum graph.
* **seamm-mdi 2026.7.15** -- consumer side: ``MDIEngine.supports(command)`` and
  ``MDIEngine.hessian()`` (returns the 3N x 3N Hessian; raises
  ``NotImplementedError`` when ``<HESSIAN`` is not advertised).
* **normal_mode_sampling_step 2026.7.15** -- the driver: gates on
  ``supports("<HESSIAN")`` and uses the analytic Hessian, else FD of the forces.

Back-pointers
-------------

* Science spec: ``~/Sites/mlff-training/2026-07-15_water-mlff-training-plan/``
  (decision #2), which now records that the sampling step and this Hessian path
  are released.
* Consumer design: the ``normal_mode_sampling_step`` campaign
  (``docs/developer_guide/campaigns/2026-07-15/`` in that repo).
