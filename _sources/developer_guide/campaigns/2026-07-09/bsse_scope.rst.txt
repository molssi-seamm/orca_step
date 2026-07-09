BSSE sub-step -- scope and design
=================================

Goal
----

A new **BSSE** sub-step of the ORCA plug-in that computes the
counterpoise-corrected (Boys--Bernardi) energy **and gradient** of a two-fragment
complex in a single ORCA run, and stores them exactly like the Energy sub-step's
``energy`` and ``gradients`` results. The immediate driver is MLFF training data
that is free of basis-set superposition error (BSSE) on both the energy surface
and the forces.

The heavy lifting is done by the ORCA *Compound* script ``BSSEGradient.cmp``
(D. G. Liakos & F. Neese, May 2024, updated June 2025). SEAMM's job is only to
(1) prepare the ghost-flagged geometry, (2) inject the method/basis/charge/mult
and options into the Compound block, (3) run ORCA once, and (4) parse the single
``result.engrad`` it writes.

What the Compound script computes
---------------------------------

For a complex partitioned into fragments A and B, it runs five calculations and
assembles the counterpoise-corrected **total** energy and gradient (not merely an
interaction energy):

.. math::

   E_\text{CP}  &= E_{AB}(AB) - \big[E_A(AB) - E_A(A)\big] - \big[E_B(AB) - E_B(B)\big] \\
   \nabla E_\text{CP} &= \nabla E_{AB}(AB) - \big[\nabla E_A(AB) - \nabla E_A(A)\big]
                       - \big[\nabla E_B(AB) - \nabla E_B(B)\big]

where ``X(AB)`` means "fragment X evaluated in the *full dimer basis*" -- i.e. the
other fragment is present as **ghost centres** (basis functions, no nucleus, no
electrons). The five calculations are: the dimer ``AB``; fragment ``A`` with B as
ghosts; ``A`` alone; fragment ``B`` with A as ghosts; ``B`` alone.

The ghost centres carry Pulay forces, so :math:`\nabla E_A(AB)` has non-zero
components on B's atoms; the script's per-atom bookkeeping (its ``CreateBSSE`` +
the ghost loop) exists precisely to map those back onto the full-complex atom
list. The final ``result.engrad`` is a standard ORCA EnGrad file over **all**
atoms of the complex, holding :math:`E_\text{CP}` and :math:`\nabla E_\text{CP}`.

How it runs today (the raw script)
----------------------------------

* Input: an xyz file (``molecule``, default ``01.xyz``) containing the **full
  complex with fragment B pre-marked as ghost atoms**. ``fragA.CreateBSSE()``
  reads the ghost flags and generates the remaining four geometries in bohr.
* Adjustable ``Variable``\ s: ``method``, ``basis``, ``restOfInput``, ``charge``,
  ``mult``, ``DoOptimization`` (relax the *free* monomers -- gives the CP
  correction relative to relaxed monomers), and file/cleanup switches.
* Output: prints all gradients and, with ``produceEnGradFile = true`` (the
  default), writes ``result.engrad``.

Design in SEAMM
---------------

Follow the standard per-capability quartet (see ``mopac_step`` / the existing
``energy`` sub-step). The BSSE node subclasses **Energy** so it inherits the
model-chemistry, method, basis, grid, SCF-convergence, and SThresh controls
unchanged.

======================  ====================================================
File                    Role
======================  ====================================================
``bsse_step.py``        Factory (entry-point target); ``my_description`` with
                        name "BSSE", group "Calculations".
``bsse.py``             ``BSSE(Energy)`` node -- overrides the run path to emit
                        a ``%Compound`` block and parse ``result.engrad``.
``bsse_parameters.py``  ``BSSEParameters(EnergyParameters)`` -- adds the
                        fragment-definition and BSSE-specific controls.
``tk_bsse.py``          ``TkBSSE(TkEnergy)`` -- adds the fragment controls to
                        the dialog.
======================  ====================================================

Register under both entry-point groups so it appears in the ORCA sub-flowchart:
``org.molssi.seamm.orca_step`` (headless) and ``org.molssi.seamm.orca_step.tk``
(GUI), alongside ``Energy`` and ``Optimization``.

Input generation (new path in ``orca_base``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The normal ``run_orca`` builds a single ``! keywords`` + inline geometry input;
the Compound job is structurally different (a ``%Compound "...end"`` block that
reads an external ``*xyzfile``). Add a sibling method, e.g.
``run_orca_compound(compound_text, xyz_files, engrad_name="result.engrad")``,
that writes the ``%pal``/``%maxcore`` preamble, the Compound block, and the
ghost-flagged xyz file(s), runs ORCA, and returns the parsed EnGrad. ``bsse.py``
builds ``compound_text`` from the ``.cmp`` template (shipped in ``data/``) by
substituting the ORCA ``with``-style variables from the resolved parameters
(reusing ``Energy.keyword_line`` for ``method``/``basis``/``restOfInput``).

Ghost geometry
~~~~~~~~~~~~~~~

The one genuinely new piece of chemistry input is writing the complex with
fragment B as ghosts. In an ORCA coordinate block a ghost centre is the element
symbol with a trailing colon, e.g. ``O:``  (**verify exact syntax against the
target ORCA version during implementation** -- the design should not hard-code
it in more than one place). ``geometry_block`` gains an optional
``ghost_atoms=<set of indices>`` argument that appends the colon for those atoms;
the BSSE node writes the complex once with B's atoms flagged.

Fragment definition (the main new UI question)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The node must know which atoms are fragment A vs B. Proposed control, defaulting
to zero-configuration for the common case:

* ``"fragments"`` = ``"auto (two molecules)"`` (default) | ``"by atom
  selection"``. Auto uses ``configuration.find_molecules()`` and requires
  exactly two molecules (error otherwise, with a clear message). "By atom
  selection" exposes an atom-index / subset selector for fragment A; the
  remainder is fragment B.
* Upstream fit: ``dimer_builder_step`` already produces two-molecule
  configurations, so ``auto`` will "just work" for the MLFF dimer campaign.

Options to expose
~~~~~~~~~~~~~~~~~

* ``DoOptimization`` (default off): relax the free monomers before the
  correction. Off is correct for a rigid-geometry PES/MLFF target; on gives the
  CP correction relative to relaxed monomers.
* Everything else (method, basis, ``restOfInput`` = ORCA "extra keywords",
  grid, SCF convergence, SThresh) is inherited from Energy.

Results
~~~~~~~

Take the corrected :math:`\nabla E_\text{CP}` (E_h/bohr) from ``result.engrad``,
but compute :math:`E_\text{CP}` (E_h) in Python from each sub-job's
``FINAL SINGLE POINT ENERGY`` rather than from ``result.engrad`` -- the script's
EnGrad energy uses ORCA's ``SCF_Energy``, which omits the MP2 correlation of a
double hybrid (see ``NOTES_implementation.rst``). Store them through the **same**
metadata result names the Energy sub-step uses (``energy`` and ``gradients``),
plus ``uncorrected energy`` and ``bsse correction``, so downstream MLFF tooling
and the property database treat a BSSE run like any other energy+gradient run.
Tag stored properties with the level of theory so corrected data is not silently
mixed with uncorrected data.

Known limitations (scope boundaries of Phase 1)
-----------------------------------------------

* **ORCA only, two fragments only.** ``CreateBSSE`` is hard-wired to two
  fragments; N-body needs the general step.
* **One charge / multiplicity for all sub-calculations.** The script passes the
  complex's ``charge``/``mult`` to the monomer steps too, so Phase 1 is valid
  only when both monomers share the complex's charge and spin (the usual
  neutral-singlet-monomer case). Charged or open-shell fragments need
  per-fragment charge/mult -- deferred to the general step. **Guard this in the
  node** and refuse other cases with a clear error rather than returning wrong
  numbers.
* **Analytic-gradient methods only** (the script always requests ``EnGrad``);
  reuse ``Energy``'s existing analytic/numeric gradient check to gate the method.
* **No MDI.** Each configuration is one ORCA Compound job (five internal SCFs);
  cost per config is ~5x a plain single point. Fine for offline training-set
  generation driven by a SEAMM Loop.

Validation
----------

* **Reproduce the raw script.** Run ``BSSEGradient.cmp`` by hand on a water dimer
  (BP86 or the campaign method) and confirm the SEAMM sub-step returns the same
  ``result.engrad`` energy and gradient to full precision.
* **Physical check.** Corrected interaction energy is less negative than the
  uncorrected one; the correction shrinks as the basis grows (e.g. def2-SVP ->
  def2-TZVPP).
* **Gradient consistency.** Finite-difference :math:`E_\text{CP}` vs the returned
  :math:`\nabla E_\text{CP}` on a few displacements.
* Add a ``tests/`` case with a small, fast dimer.

Citations
---------

Add to the step's references: Boys & Bernardi, *Mol. Phys.* **19**, 553 (1970)
for the counterpoise method; the ORCA *Compound* feature and the
``BSSEGradient.cmp`` authors (Liakos & Neese). Keep the existing ORCA / method /
basis citations.

Implementation task list
------------------------

#. Ship ``BSSEGradient.cmp`` as a template in ``orca_step/data/`` (parameterized
   via ORCA ``with`` variables).
#. ``orca_base``: ``geometry_block(..., ghost_atoms=None)`` and
   ``run_orca_compound(...)`` + EnGrad parsing reuse.
#. ``bsse_parameters.py`` (fragments, DoOptimization) and ``tk_bsse.py``.
#. ``bsse.py``: resolve fragments, build the Compound text, run, parse, store
   ``energy``/``gradients``; guard the charge/mult and gradient-availability
   limits.
#. Entry points in ``setup.py`` (headless + tk) and re-exports in
   ``__init__.py``.
#. Citations, ``tests/``, user-guide section, HISTORY entry; release via the
   ``release-seamm-plugin`` skill.
