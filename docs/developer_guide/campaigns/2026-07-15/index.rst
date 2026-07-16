2026-07-15 Frequencies sub-step and the <HESSIAN MDI command
============================================================

Add an ORCA **Frequencies** sub-step -- the Hessian, harmonic vibrational
frequencies, IR intensities, and thermochemistry via ORCA's analytic
(``AnFreq``) or numerical (``NumFreq``) second derivatives -- and a custom
**``<HESSIAN``** MDI command that serves the analytic Cartesian Hessian over a
warm MDI connection.

This is the ORCA end of a cross-repo *Hessian-over-MDI* effort for the MLFF
reference-data campaign: a driver that needs a Hessian gets an analytic one
straight from the resident engine when the method supports it, and
finite-differences the forces otherwise. The first consumer is the
``normal_mode_sampling_step`` plug-in, which reaches the engine through
``seamm_mdi.MDIEngine``. The **science spec** lives in the ``~/Sites`` lab
notebook (the water |rarr| electrolyte MLFF training-set plan, decision #2,
"Internal sampling = Wigner normal-mode sampling"); this page is the ORCA
implementation record.

Released in orca_step **2026.7.15** (the feature) with follow-up polish in
**2026.7.16**.

.. |rarr| unicode:: U+2192

Contents:

.. toctree::
   :glob:
   :maxdepth: 2

   *scope*
   NOTES_*
