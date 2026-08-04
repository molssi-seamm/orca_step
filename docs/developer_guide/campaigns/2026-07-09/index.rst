2026-07-09 BSSE (counterpoise) sub-step
=======================================

Add a **BSSE** sub-step to the ORCA plug-in that returns the
counterpoise-corrected energy and gradient of a molecular complex, for use as
machine-learned-force-field (MLFF) training data. It wraps the ORCA *Compound*
script ``BSSEGradient.cmp`` (D. G. Liakos & F. Neese, 2024/2025).

This is the fast, ORCA-specific first step (two fragments, neutral singlet
only, driven by a single ORCA *Compound* process). It was later generalized
to an arbitrary number of fragments with independent per-fragment charge,
backed by a new shared library rather than a second ORCA-specific special
case -- see :doc:`the 2026-08-03 campaign <../2026-08-03/index>`, whose
canonical design document lives in the sibling ``seamm_bsse`` repository.
This Compound-script path is kept as that generalization's permanent N = 2
regression reference, not deleted, even though ``run()`` no longer calls it.

Contents:

.. toctree::
   :glob:
   :maxdepth: 2

   *scope*
   NOTES_*
