2026-07-09 BSSE (counterpoise) sub-step
=======================================

Add a **BSSE** sub-step to the ORCA plug-in that returns the
counterpoise-corrected energy and gradient of a molecular complex, for use as
machine-learned-force-field (MLFF) training data. It wraps the ORCA *Compound*
script ``BSSEGradient.cmp`` (D. G. Liakos & F. Neese, 2024/2025).

This is the fast, ORCA-specific first step. The code-agnostic, N-fragment
generalization is a separate, larger project with its own design document (see
the general ``bsse_step`` plug-in).

Contents:

.. toctree::
   :glob:
   :maxdepth: 2

   *scope*
   NOTES_*
