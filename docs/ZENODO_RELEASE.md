# Zenodo and GitHub release checklist

1. Resolve the two implementation choices in `REPRODUCIBILITY.md` and update the manuscript/repository consistently.
2. Confirm MIT (code) and CC BY 4.0 (data/figures) with every author.
3. Add missing author ORCIDs and exact affiliations to `.zenodo.json` and `CITATION.cff`.
4. Run `python -m pytest` and reproduce all available figures.
5. Rebuild `data/MANIFEST.csv` and `checksums.sha256` after any data change.
6. Tag the GitHub release, preferably `v1.0.0` for the accepted-paper archive.
7. Enable the GitHub–Zenodo integration and archive that release.
8. Replace the temporary Zenodo badge and complete the citation metadata with the concept DOI and version DOI.
9. After ApJ publication, add journal DOI, volume, article number, ADS bibcode, and the bidirectional related identifier in Zenodo.
10. Cite the immutable Zenodo DOI in the paper’s final software/data statement.
