#!/usr/bin/env python3
"""Merge array-job result tables while preserving one commented header."""
from __future__ import annotations

import argparse
from pathlib import Path

from moon_packing.io import read_results, write_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = []
    n_moons = None
    for path in args.inputs:
        table = read_results(path)
        n_moons = table.n_moons if n_moons is None else n_moons
        if table.n_moons != n_moons:
            raise ValueError("All tiles must have the same number of moons")
        for j in range(len(table)):
            row = {name: table[name][j] for name in table.numeric}
            row["fate"] = table.fate[j]
            rows.append(row)
    if n_moons is None:
        raise ValueError("No input data")
    write_results(args.output, rows, n_moons)
    print(args.output)


if __name__ == "__main__":
    main()
