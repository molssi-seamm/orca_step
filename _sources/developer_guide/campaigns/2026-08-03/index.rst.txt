2026-08-03 N-fragment, charge-aware counterpoise (seamm_bsse)
================================================================

Generalizes the 2026-07-09 **BSSE** sub-step (two fragments, neutral singlet
only, one ORCA *Compound* process) to an arbitrary number of fragments with
independent per-fragment charge, backed by a new shared library,
**seamm_bsse**, rather than a second ORCA-specific special case.

The canonical design document -- the physics, the architecture (why a
separate library plus a thin per-engine sub-step rather than one
code-agnostic driver step), and the confirmed ORCA ghost-atom indexing
property that makes the library's ``combine()`` padding-only -- lives in the
sibling ``seamm_bsse`` repository:
``seamm_bsse/docs/developer_guide/campaigns/2026-08-03/bsse_architecture.rst``.
This page is the ORCA-side implementation record and validation log.

What changed here
------------------

* ``orca_base.py``: ``run_orca_job`` extracted from ``run_orca`` --
  geometry/charge/multiplicity/directory now come from arguments instead of
  ``self.get_system_configuration()``/``self.directory``, since a
  counterpoise sub-job differs from the node's own system in exactly those
  ways. ``run_orca`` is now a thin wrapper. ``geometry_block`` gained
  ``atom_indices=``/``ghost_atoms=`` (backward compatible).
* ``bsse.py``: ``_fragments()`` (renamed from ``_fragment_atoms``) returns a
  list of ``seamm_bsse.Fragment`` for N groups with independent charge;
  ``run()`` generates the 2N + 1 job specs via
  ``seamm_bsse.generate_job_specs``, runs each through ``run_orca_job`` in
  its own sub-directory, and assembles the result via
  ``seamm_bsse.combine()``.
* ``bsse_parameters.py``/``tk_bsse.py``: ``"fragments"`` gained
  ``"auto (molecules)"`` (any number of molecules, not just two);
  ``"fragment A atoms"`` became ``"fragment atoms"`` (semicolon-separated
  groups, one per fragment); new ``"fragment charges"`` parameter.
* The original ORCA *Compound* path (``bssegradient.cmp``/``bssenergy.cmp``,
  ``_compound_input``/``run_orca_compound``/``_parse_compound_energies`` in
  ``bsse.py``) is **not removed** -- it stays as the permanent N = 2
  regression reference the new path was validated against, just no longer
  called from ``run()``.

Validation
----------

Three real-ORCA (not stubbed) validation scripts, each building a real
``molsystem`` Configuration and a minimal-but-real execution harness
(``seamm_exec``'s local executor, actual ORCA 6.1.1) rather than a
``.flow``/``run_flowchart`` -- needed because the old Compound path is no
longer reachable through a live flowchart run, so both paths (old, called
directly via the still-present helper methods; new, plain ``node.run()``)
have to run side by side in the same script to be compared.

:download:`validate_bsse_m2.py <validate_bsse_m2.py>`
   N = 2, neutral-singlet regression gate: does the new path reproduce the
   old Compound-script path? **Passed** on the water dimer
   (``H2O-H2O.sdf``), FEC-FEC (``FEC-FEC dimer opt PM6-ORG.sdf``, incl. F),
   and EC-EC (first record of the production ``EC_dimers_run_1.sdf``, an
   NMS-displaced sample, not a hand-picked equilibrium geometry) at
   HF/def2-SVP. Energies agreed to ~2-7e-9 E\ :sub:`h`, gradients to
   ~5.6-6.9e-8 E\ :sub:`h`/bohr on all three -- SCF-noise level.

:download:`validate_bsse_m3.py <validate_bsse_m3.py>`
   Per-fragment charge on real charged two-body systems, B3LYP-D3BJ/def2-TZVP,
   vs. approximate literature binding energies. Na\ :sup:`+`\ ···Cl\
   :sup:`-` well minimum &minus;136 kcal/mol at R &asymp; 2.4-2.6 Å (lit.
   &asymp; &minus;133 at R\ :sub:`e` &asymp; 2.36 Å, from a Born-Haber cycle
   over atomic D0/IP/EA), decaying smoothly to &minus;64 kcal/mol by 7 Å
   (points from the already-built 40-point ``Na-Cl.sdf`` R-scan); Na\
   :sup:`+`\ ···H\ :sub:`2`\ O = &minus;26.1 kcal/mol (lit. &asymp;
   &minus;24); Cl\ :sup:`-`\ ···H\ :sub:`2`\ O = &minus;15.6 kcal/mol (lit.
   &asymp; &minus;13). All within a few kcal/mol on unoptimized,
   literature-informed geometries -- the expected residual from not
   relaxing geometry, not a wiring defect.

:download:`validate_bsse_n3.py <validate_bsse_n3.py>`
   N = 3 on real ORCA: a hand-built Na\ :sup:`+`\ ···Cl\ :sup:`-`\ ···H\
   :sub:`2`\ O trimer (contact ion pair plus a water coordinating Na\
   :sup:`+`\  from a different direction than Cl\ :sup:`-`), same
   B3LYP-D3BJ/def2-TZVP level. Ran the full 2N + 1 = 7 jobs correctly
   (per-fragment charges +1/&minus;1/0 confirmed in the printed output). CP
   interaction &minus;150.2 kcal/mol: deeper than the Na\ :sup:`+`\ ···Cl\
   :sup:`-` pair alone (&minus;136) but *less* than the naive pairwise sum
   (&minus;136 + &minus;26 &asymp; &minus;162) -- the expected
   cooperative-saturation non-additivity, not a bug.

:download:`optimize_ion_water.py <optimize_ion_water.py>`
   Geometry optimization for the Na\ :sup:`+`\ ···H\ :sub:`2`\ O and
   Cl\ :sup:`-`\ ···H\ :sub:`2`\ O pairs (B3LYP-D3BJ/def2-TZVP, ``TightOpt``,
   from M3's literature-informed starting geometries -- via ``orca_step``'s
   own ``Optimization`` sub-step, then ``BSSE`` at the relaxed geometry).
   Na-O tightens from the 2.30 Å starting guess to ~2.21 Å. CP interaction
   energy at the true minimum: Na\ :sup:`+`\ ···H\ :sub:`2`\ O = &minus;26.4
   kcal/mol (M3 unoptimized: &minus;26.1; lit. &asymp; &minus;24); Cl\
   :sup:`-`\ ···H\ :sub:`2`\ O = &minus;17.5 kcal/mol (M3 unoptimized:
   &minus;15.6; lit. &asymp; &minus;13). Optimizing moved *both* numbers
   further from the approximate literature references, not closer --
   traced to a basis-set effect, not a geometry or CP-wiring problem (see
   below).

:download:`check_diffuse_basis.py <check_diffuse_basis.py>`
   Same optimized geometries, ``def2-TZVP`` vs. ``def2-TZVPPD`` (the
   campaign's actual production basis, with diffuse functions) as a
   single-point ``BSSE`` re-evaluation -- isolates the basis effect from any
   geometry-relaxation effect. Diffuse functions close roughly half the gap
   to literature in both cases: Cl\ :sup:`-`\ ···H\ :sub:`2`\ O &minus;17.5
   -> &minus;15.5 kcal/mol (lit. &asymp; &minus;13); Na\ :sup:`+`\ ···H\
   :sub:`2`\ O &minus;26.4 -> &minus;25.1 kcal/mol (lit. &asymp; &minus;24).
   Confirms the residual is the same "diffuse functions matter for
   dispersion/diffuse-anion interactions" finding the campaign already made
   for neutral dimers (``def2-TZVPPD`` vs. plain ``def2-TZVPP``,
   see ``bsse-corrected-gradients`` era notes) -- ``def2-TZVP`` (no diffuse
   "D") understates ion-water binding, most visibly for the more diffuse
   Cl\ :sup:`-`\ anion. The remaining ~1-2.5 kcal/mol gap after adding
   diffuse functions is plausible from (a) the approximate literature
   numbers being ballpark, not rigorous CCSD(T)/CBS benchmarks, (b)
   B3LYP-D3BJ vs. the campaign's actual production method
   (revDSD-PBEP86-D4), and (c) no re-optimization at the diffuse basis
   (only a single point at the ``def2-TZVP``-optimized geometry).

Not done
--------

* The Psi4 sub-step (the cross-engine check N = 3 was originally paired
  with) is done -- see the sibling ``psi4_step`` repository's own
  ``docs/developer_guide/campaigns/2026-08-04/``.
* A systematic angular scan for the Na\ :sup:`+`\ ···H\ :sub:`2`\ O and
  Cl\ :sup:`-`\ ···H\ :sub:`2`\ O pairs (the S66/angular-test pattern) --
  only the single equilibrium geometry has been optimized so far.
* Re-optimizing at ``def2-TZVPPD`` (or the actual production method,
  revDSD-PBEP86-D4) rather than reusing the ``def2-TZVP``-optimized
  geometry as a single-point re-evaluation.
