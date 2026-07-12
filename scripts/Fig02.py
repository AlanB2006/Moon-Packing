#!/usr/bin/env python3

from pathlib import Path

from moon_packing.plotting import plot_eccentricity_map

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "figures" / "reproduced" / "figure08_luna_eccentricity_map.pdf"


def main() -> int:
    output = plot_eccentricity_map("Luna", root=ROOT, output=OUTPUT)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
