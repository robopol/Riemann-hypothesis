from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "required_prime_harmonic_upper_rows.csv"
DEFAULT_JSON = ROOT / "required_harmonic_constant_report.json"
DEFAULT_CSV = ROOT / "required_harmonic_constant_rows.csv"

MEISSEL_MERTENS_PRIME_CONSTANT = 0.261497212847642783755426838608695859


@dataclass
class RequiredHarmonicConstantRow:
    x_prime: int
    actual_constant: float
    required_constant: float
    constant_headroom: float
    rosser_constant: float
    rosser_excess: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scale the required A(x) upper envelope as B+C/(sqrt(x)log x)."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        input_rows = list(csv.DictReader(handle))
    if not input_rows:
        raise ValueError(f"No rows found in {args.input}")

    rows: list[RequiredHarmonicConstantRow] = []
    for raw in input_rows:
        x = int(raw["x_prime"])
        scale = math.sqrt(x) * math.log(x)
        actual_constant = (float(raw["actual_a_x"]) - MEISSEL_MERTENS_PRIME_CONSTANT) * scale
        required_constant = (
            float(raw["required_a_upper"]) - MEISSEL_MERTENS_PRIME_CONSTANT
        ) * scale
        rosser_constant = 2.0
        rows.append(
            RequiredHarmonicConstantRow(
                x_prime=x,
                actual_constant=actual_constant,
                required_constant=required_constant,
                constant_headroom=required_constant - actual_constant,
                rosser_constant=rosser_constant,
                rosser_excess=rosser_constant - required_constant,
            )
        )

    payload = {
        "notes": [
            "This rescales the required A(x) upper envelope as A(x) <= B + C/(sqrt(x) log x).",
            "The Rosser-Schoenfeld finite-range diagnostic corresponds to C=2.",
            "Positive rosser_excess means C=2 is too loose compared with the step requirement.",
            "Positive constant_headroom is the real numerical gap between exact A(x) and the required envelope.",
        ],
        "max_actual_constant": max(row.actual_constant for row in rows),
        "max_required_constant": max(row.required_constant for row in rows),
        "min_required_constant": min(row.required_constant for row in rows),
        "min_constant_headroom": min(row.constant_headroom for row in rows),
        "min_rosser_excess": min(row.rosser_excess for row in rows),
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))

    print(f"read {args.input}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Required C in A(x) <= B + C/(sqrt(x) log x):")
    for row in rows:
        print(
            f"  x={row.x_prime:>9} "
            f"C_actual={row.actual_constant:.6f} "
            f"C_required={row.required_constant:.6f} "
            f"headroom={row.constant_headroom:+.6f} "
            f"C2_excess={row.rosser_excess:+.6f}"
        )


if __name__ == "__main__":
    main()
