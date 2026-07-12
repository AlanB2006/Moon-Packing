from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass
class ResultsTable:
    path: Path
    names: list[str]
    numeric: dict[str, np.ndarray]
    fate: np.ndarray
    n_moons: int

    def __len__(self) -> int:
        return len(self.fate)

    def __getitem__(self, name: str) -> np.ndarray:
        if name == "fate":
            return self.fate
        return self.numeric[name]

    def eccentricity_range(self, moon_index: int) -> np.ndarray:
        return self[f"eM{moon_index}_max"] - self[f"eM{moon_index}_min"]


def _clean_header(line: str) -> list[str]:
    line = line.lstrip("#").strip()
    return [piece.strip() for piece in line.split(",") if piece.strip()]


def read_results(path: str | Path) -> ResultsTable:
    """Read a published result table, tolerating whitespace and mixed precision."""
    path = Path(path)
    with path.open(encoding="utf-8") as handle:
        first = handle.readline()
        if not first.startswith("#"):
            raise ValueError(f"{path} does not begin with the expected commented header")
        names = _clean_header(first)
        rows: list[list[str]] = []
        reader = csv.reader(handle, skipinitialspace=True)
        for row in reader:
            cleaned = [cell.strip() for cell in row]
            if cleaned and any(cleaned):
                rows.append(cleaned)
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    width = len(names)
    bad = [i + 2 for i, row in enumerate(rows) if len(row) != width]
    if bad:
        raise ValueError(f"Malformed rows in {path}: {bad[:5]}")
    numeric: dict[str, np.ndarray] = {}
    fate_index = names.index("fate")
    for j, name in enumerate(names):
        if j == fate_index:
            continue
        numeric[name] = np.asarray([float(row[j]) for row in rows], dtype=float)
    fate = np.asarray([row[fate_index] for row in rows], dtype=str)
    n_moons = sum(name.endswith("_max") and name.startswith("eM") for name in names)
    return ResultsTable(path, names, numeric, fate, n_moons)


def result_header(n_moons: int) -> list[str]:
    names = ["dH", "e0"]
    for i in range(1, n_moons + 1):
        names.extend([f"eM{i}_max", f"eM{i}_min", f"dA{i}_over_A0"])
    names.extend(["t_end(PM)", "wall_time", "fate"])
    return names


def write_results(path: str | Path, rows: Iterable[dict[str, float | str]], n_moons: int) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = result_header(n_moons)
    rows = list(rows)
    rows.sort(key=lambda r: (float(r["e0"]), float(r["dH"])))
    with path.open("w", newline="", encoding="utf-8") as handle:
        handle.write("# " + ", ".join(names) + "\n")
        writer = csv.writer(handle)
        for row in rows:
            out: list[str] = []
            for name in names:
                value = row[name]
                if name == "fate":
                    out.append(str(value))
                elif name in {"dH", "e0"}:
                    out.append(f"{float(value):.6g}")
                elif name == "wall_time":
                    out.append(f"{float(value):.5f}")
                else:
                    out.append(f"{float(value):.8e}")
            writer.writerow(out)
