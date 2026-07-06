=======
History
=======

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
    * The library-path directory is now exported inside the run command so it
      survives macOS System Integrity Protection (which strips DYLD_* variables
      passed through the shell); without this, parallel runs on macOS could not
      find ORCA's OpenMPI library (libmpi).
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
