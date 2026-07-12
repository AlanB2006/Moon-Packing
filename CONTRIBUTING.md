# Contributing

For scientific changes, open an issue describing the numerical or plotting
change and whether it affects archived results. Keep the archived files in
`data/`, `figures/reference/`, and `legacy/scripts/` unchanged. Add tests for
new reusable code and run:

```bash
python -m pytest
moon-packing-figures --figures 2 3 4 8 9 10
```

Changes that alter initial conditions, stopping criteria, units, dependency
versions, or table schemas must also update `docs/REPRODUCIBILITY.md` and the
Zenodo release notes.
