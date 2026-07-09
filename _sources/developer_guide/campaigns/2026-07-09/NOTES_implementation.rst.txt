Implementation notes
=====================

2026-07-09 -- built and validated
---------------------------------

Implemented as ``bsse.py`` / ``bsse_parameters.py`` / ``tk_bsse.py`` /
``bsse_step.py`` (+ entry points, ``data/bssegradient.cmp``). ``orca_base`` grew
``_resources``/``_mpi_env`` (extracted from ``run_orca``), ``run_orca_compound``,
and ``_read_engrad``.

Double-hybrid energy gotcha (important, and relevant to the general step)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The shipped Compound script reads ORCA's ``SCF_Energy`` property for each
sub-calculation. For a **double hybrid** (the campaign's ``REVDSD-PBEP86-D4/2021``)
that is only the SCF part -- it omits the (scaled) MP2 correlation, which for
Ar\ :sub:`2`/def2-TZVP is ~0.2 E\ :sub:`h`. So the script's ``result.engrad``
*energy* is an SCF-level counterpoise energy, and comparing it to the true total
gave a nonsensical "correction" of 130 kcal/mol on the first run.

The script's **gradient** is fine: it operates on ``Nuclear_Gradient``
(``Property_Base=true``), which is the full-method (relaxed-density MP2 +
dispersion) gradient.

Fix (method-general): **ignore the script's energy; compute the corrected energy
in Python** from each sub-job's ``FINAL SINGLE POINT ENERGY`` -- the true total
that matches the gradient. The five sub-jobs are delimited in ``orca.out`` by
``COMPOUND JOB 1..5`` (in the order fragA(AB), monA, fragB(AB), monB, dimer);
take the *last* ``FINAL SINGLE POINT ENERGY`` in each block so an optimized
monomer's converged value wins. Dispersion has no BSSE (ghost centres have no
nuclei) and cancels in the correction terms, so ``FINAL SINGLE POINT ENERGY`` is
consistent for ``-D`` methods too. This removed the need for an HF/DFT-only
restriction: **any analytic-gradient method works** (HF, DFT incl.
dispersion-corrected and double hybrids, MP2/RI-MP2); only numerical-gradient
methods ((DLPNO-)CCSD(T)) are refused.

Validation
~~~~~~~~~~~

``Testing/bsse.flow`` on Ar-Ar (REVDSD-PBEP86-D4/2021/def2-TZVP) gives:

* uncorrected -1054.38210538 E_h (matches a plain single-point Energy run),
* corrected   -1054.38198513 E_h,
* BSSE correction 0.00012024 E_h = **0.0755 kcal/mol** (physically sane for a
  rare-gas dimer at a triple-zeta double hybrid).

Confirmed on this machine that the ``%Compound "file.cmp" with <var=val;> end``
injection syntax and the ghost notation (element + ``:`` in the ``*xyzfile``)
both work as written.

Still open / deferred to the general step: per-fragment charge & multiplicity
(only neutral-singlet now), N > 2 fragments, non-ORCA engines, and a
finite-difference check of the corrected gradient.
