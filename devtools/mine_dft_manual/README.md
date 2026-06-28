# Mining DFT functional citations from the ORCA manual

`parse_dft_manual.py` builds the DFT-functional citation data shipped in the
ORCA step:

- appends bibtex entries to `orca_step/data/references.bib` (keys `orca_dft_<id>`)
- writes `orca_step/data/dft_functionals.json` = `{functional keyword: [keys]}`
  (keys only — the bibtex lives only in references.bib)

## To regenerate (e.g. when ORCA updates the manual)

Download the two source pages into this directory, then run the script:

```bash
curl -sL -o dft.html \
  https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/DensityFunctionalTheory.html
curl -sL -o manual_index.html \
  https://orca-manual.mpi-muelheim.mpg.de/index.html
python parse_dft_manual.py
```

The DFT page provides `functional keyword -> reference number`; the manual
`index.html` is the global bibliography (`reference id -> citation text + DOI`).
Functional tables use the header "Keyword" (hybrids) or "Keywords"
(local/double-hybrid). Some older references have no DOI (kept as citation text).
