=======
History
=======

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
