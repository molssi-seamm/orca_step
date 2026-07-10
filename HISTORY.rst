=======
History
=======

2026.7.10 -- BSSE: energy-only mode (enables CCSD(T))
    * New **Compute the gradient** control on the BSSE sub-step. With it set to
      ``no`` (energy only), the counterpoise correction runs without a gradient,
      which is cheaper and works for methods that have no analytic gradient in
      ORCA -- notably ``CCSD(T)`` / ``DLPNO-CCSD(T)`` -- for gold-standard
      counterpoise interaction energies. The default (``yes``) is unchanged and
      still produces the corrected energy and gradient for MLFF training.

2026.7.9.2 -- BSSE (counterpoise) sub-step
    * New **BSSE** sub-step: the counterpoise-corrected (Boys--Bernardi) energy
      and gradient of a two-fragment complex, in a single ORCA run, for
      BSSE-free machine-learned-force-field training data. It drives ORCA's
      Compound facility (the BSSEGradient script by D. G. Liakos & F. Neese).
    * Reports the BSSE-corrected energy, the uncorrected (raw) energy, and the
      correction (in E_h and kcal/mol), plus the corrected gradient; each can be
      saved from the Results tab.
    * Fragments are found automatically from the two molecules in the structure
      (so it works directly on a Dimer Builder dimer) or specified by atom; an
      option relaxes the free monomers before the correction.
    * First version: a neutral, closed-shell complex of exactly two fragments
      with an ORCA-internal basis set. Any analytic-gradient method works,
      including dispersion-corrected and double-hybrid DFT and MP2.

2026.7.9.1 -- SCF SThresh control on the Energy step
    * New **SCF SThresh** control: set ORCA's SCF convergence threshold (the
      ``%scf SThresh`` value, in E_h). Leave it at ``default`` to let the SCF
      convergence preset (or ORCA's own default) govern SThresh, or give an
      explicit value -- e.g. ORCA's nominal ``1.0e-07`` -- to write it out and
      override whatever the preset would otherwise set. Lower it for a tighter
      SCF (smoother energies and forces), raise it to converge more loosely.
      Available on the Energy and Optimization sub-steps.

2026.7.9 -- MDI engine, integration-grid and SCF controls, and config fixes
    * ORCA can now be driven as a persistent MDI engine, so steps that set up a
      model chemistry and evaluate it at many geometries (for example the Dimer
      Builder's energy-based contact search) can use ORCA -- for methods with an
      analytic gradient. Set this up with a Model Chemistry step; the ORCA step
      does not need configuring for it.
    * New **Integration grid** control: choose ORCA's grid preset (DEFGRID1,
      DEFGRID2, or DEFGRID3), or leave ORCA's default.
    * New **SCF convergence** control: choose the convergence-tolerance preset
      (SLOPPYSCF ... EXTREMESCF), or leave ORCA's default. It defaults to
      TIGHTSCF for smooth energies and forces (previously this was applied via
      the extra-keywords default, which is now empty).
    * Bugfix: a ``$variable`` typed into the basis-set field (e.g. to vary the
      basis in a Loop) is now expanded to its value instead of being passed to
      ORCA literally.
    * Bugfix: the pre-run description of a step showed the basis as a raw
      dictionary; it now shows the basis name.
    * How to find and launch ORCA -- its executable path and, for parallel runs,
      the OpenMPI library directory -- now lives in ``~/SEAMM/orca.ini`` (a
      template is created on install). The ``[orca-step]`` section of the main
      SEAMM configuration keeps only the user run options (ncores, memory).
      NOTE: if you previously set ``library-path`` (or a path) in the
      ``[orca-step]`` section, move it into ``~/SEAMM/orca.ini``.

2026.7.8 -- Ordered basis-set list and complete-basis-set (CBS) extrapolation
    * The basis-set list is now ordered by family and, within a family, into
      valence / polarization / diffuse ladders that each rise DZ -> TZ -> QZ ->
      5Z, so a sensible progression is a single ladder read top to bottom.
    * New complete-basis-set (CBS) extrapolation on the Energy step: set
      'Basis-set extrapolation' to 2/3, 3/4, or 4/5 and pick a family (cc,
      aug-cc, def2, or ANO). This is a single ORCA job (its Extrapolate keyword)
      that runs both basis sets and extrapolates the SCF and correlation parts.
      When it is on, the fixed basis set is ignored, and gradients are not
      available (ORCA has no gradient for an extrapolated energy).
    * The CBS control is hidden for the Optimization step, which needs a
      gradient that an extrapolated energy does not provide.

2026.7.6.1 -- Bugfix: parallel execution and DFT functionals via Model Chemistry
    * Bugfix: the library-path (and orca-path) settings were read under the
      wrong key and so were ignored; they are now applied. The matching OpenMPI
      'mpirun' (the sibling 'bin' of library-path) is put on PATH so ORCA
      launches its workers with the correct OpenMPI -- a mismatched one (e.g. a
      newer system OpenMPI) causes parallel runs to abort with a BLAS-ERROR.
      The loader variables are also exported inside the run command so they
      survive macOS System Integrity Protection. (On macOS ORCA does not pass
      DYLD_* to its MPI sub-processes, so the OpenMPI libraries must additionally
      be on the default loader path, e.g. symlinked into /usr/local/lib; see the
      User Guide.)
    * Bugfix: the DFT functionals are again selectable through the Model
      Chemistry step (each functional is offered as a method); this regressed
      when the functionals moved out of the method list. Keywords containing '/'
      (e.g. REVDSD-PBEP86-D4/2021) appear with '_' in the model-chemistry string,
      since '/' is reserved there, and are translated back to the real keyword
      when the calculation runs.

2026.7.6 -- All ORCA functionals, forces, database properties, and parallel execution
    * Density functional theory now offers the complete set of ORCA functionals
      (117 of them), organized by type: pick a functional type (local, GGA,
      meta-GGA, hybrid, range-separated hybrid, or double-hybrid) and then the
      functional itself, including the double hybrids such as
      REVDSD-PBEP86-D4/2021.
    * Gradients (forces) are produced with the correct ORCA method automatically:
      the analytic gradient where ORCA has one, or the numerical gradient where
      it does not (for example DLPNO-CCSD(T), and the non-self-consistent
      wB97M(2) and wB97X-2 functionals). A note is printed when the slower
      numerical gradient is used.
    * Results can now be saved to the property database, including the gradient,
      the dipole-moment vector, the Mulliken, Löwdin, and Hirshfeld charges, the
      Mayer valences, and the rotational constants, in addition to the energies
      and other scalar results.
    * The curated basis-set list has been filled out across the Pople, Dunning
      (correlation-consistent), and Karlsruhe def2 families, including their
      diffuse and minimally-augmented variants.
    * Choosing 'Basis Set Exchange' as the basis-set source now opens the picker
      directly (it remains available from the '...' button as well).
    * ORCA runs in parallel by default, using the cores the machine or batch job
      provides. The number of cores and the memory per process can be set in the
      [orca-step] section of orca.ini (the ncores and memory options); parallel
      runs need ORCA's OpenMPI runtime, whose location can be given with
      library-path.
    * Documentation: a full User Guide covering methods and functionals, basis
      sets, forces, saving results, and parallel execution.

2026.6.28.1 -- A Basis Set Exchange basis-set picker
    * The basis set now uses the shared Basis Set Exchange picker: type a name,
      pick a common one from the list, or press '...' to browse any basis from
      the Exchange, filtered to the elements you select on a periodic table. A
      choice from the Exchange is stored as 'bse:NAME', and the element selection
      is remembered so the picker is restored when the flowchart is reopened.
    * Bugfix: the basis-set source control no longer appears on its own when the
      model chemistry is used; it is shown only with an explicit method and basis.

2026.6.28 -- Properties, gradients, citations, and wavefunction export
    * Reports many properties from a single calculation: HOMO/LUMO (and the next
      orbitals) and the gap, the dipole moment, rotational constants, <S^2>, the
      Mulliken, Löwdin, and Hirshfeld atomic charges, the Mayer bond orders and
      valences, and the optional dipole polarizability.
    * The Mayer bond orders and Hirshfeld charges can be written to a CSV file
      and applied to the structure.
    * Energy gradients can be requested and are written to Results.json for use by
      driver steps such as Thermochemistry and Reaction Path.
    * Full citations for each run: the ORCA program, the DFT functional (from the
      ORCA manual), the basis set (via the Basis Set Exchange), and the supporting
      integral and exchange-correlation libraries.
    * Basis sets can be taken from the Basis Set Exchange, including a 'bse:NAME'
      shorthand that forces a single basis from the Exchange.
    * Can write an analytic wavefunction (.wfx, via orca_2aim) for a following
      Atomic Charges step to partition into DDEC6 charges.
    * Fixed: the Results tab in the GUI was empty; it now lists the available
      results to save to variables, tables, or JSON.

2026.6.27 -- Initial release of the ORCA step
    * A sub-flowchart ORCA plug-in with Energy and Optimization sub-steps.
    * Single-point energies and geometry optimizations, including DLPNO-CCSD(T).
    * The method and basis set can be set explicitly, or taken from a preceding
      Model Chemistry step. Basis sets use ORCA's built-in families.
