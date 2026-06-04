from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "required_harmonic_constant_rows.csv"
DEFAULT_JSON = ROOT / "literature_harmonic_bounds_report.json"
DEFAULT_CSV = ROOT / "literature_harmonic_bounds_rows.csv"


@dataclass
class LiteratureHarmonicBoundRow:
    x_prime: int
    actual_constant: float
    required_constant: float
    target_constant_1_margin: float
    rosser_schoenfeld_constant: float
    rosser_schoenfeld_margin: float
    dusart_theorem_6_10_constant: float
    dusart_margin: float
    optimistic_0_2_log3_constant: float
    optimistic_0_2_log3_margin: float
    axler_2018_constant: float | None
    axler_2018_margin: float | None
    rh_reciprocal_constant_3_over_8pi: float
    rh_reciprocal_margin: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare required prime-harmonic constants with literature-scale bounds."
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

    rows: list[LiteratureHarmonicBoundRow] = []
    for raw in input_rows:
        x = int(raw["x_prime"])
        log_x = math.log(x)
        required = float(raw["required_constant"])
        actual = float(raw["actual_constant"])

        # Rosser-Schoenfeld finite-window diagnostic scale:
        # |A(x)-B| <= 2/(sqrt(x) log x), i.e. C=2.
        rosser_c = 2.0

        # Dusart Theorem 6.10:
        # A(x)-B <= 1/(10 log^2 x) + 4/(15 log^3 x), for x > 10372.
        # In the C/(sqrt(x)log(x)) scale this corresponds to
        # C = sqrt(x)log(x) * (1/(10log^2 x) + 4/(15log^3 x)).
        dusart_c = math.sqrt(x) * (
            0.1 / log_x + (4.0 / 15.0) / (log_x**2)
        )

        # Dusart's later 0.2/log^3 reciprocal-prime bound, converted to the
        # common C/(sqrt(x)log(x)) scale.
        optimistic_c = 0.2 * math.sqrt(x) / (log_x**2)

        # Axler 2018 upper bound:
        # A(x)-B <= 1/(20 log^3 x) + 3/(16 log^4 x), for x >= 46,909,074.
        axler_c = None
        axler_margin = None
        if x >= 46_909_074:
            axler_c = math.sqrt(x) * (
                1.0 / (20.0 * log_x**2) + 3.0 / (16.0 * log_x**3)
            )
            axler_margin = required - axler_c

        # RH-conditional reciprocal-prime bound reported in the
        # Schoenfeld/Lee-Nosal family:
        # |A(x)-B| <= (3 log(x)+4)/(8*pi*sqrt(x)).
        # In our scale this is C = (3 log^2(x)+4log(x))/(8*pi).
        rh_c = (3.0 * (log_x**2) + 4.0 * log_x) / (8.0 * math.pi)

        rows.append(
            LiteratureHarmonicBoundRow(
                x_prime=x,
                actual_constant=actual,
                required_constant=required,
                target_constant_1_margin=required - 1.0,
                rosser_schoenfeld_constant=rosser_c,
                rosser_schoenfeld_margin=required - rosser_c,
                dusart_theorem_6_10_constant=dusart_c,
                dusart_margin=required - dusart_c,
                optimistic_0_2_log3_constant=optimistic_c,
                optimistic_0_2_log3_margin=required - optimistic_c,
                axler_2018_constant=axler_c,
                axler_2018_margin=axler_margin,
                rh_reciprocal_constant_3_over_8pi=rh_c,
                rh_reciprocal_margin=required - rh_c,
            )
        )

    payload = {
        "notes": [
            "All bounds are converted to the common scale A(x)-B <= C/(sqrt(x)log x).",
            "Positive margin means the bound is strong enough for the required ledger constant at that endpoint.",
            "Dusart 2010 Theorem 6.10 is rigorous but much too large in this scale.",
            "The 0.2/log^3 entry is Dusart's later reciprocal-prime bound, not an analytic closure by itself.",
            "Axler 2018 gives a sharper upper bound from x >= 46,909,074 and can close some finite endpoints.",
            "The RH-conditional (3log(x)+4)/(8pi sqrt(x)) bound is also too large in this scale.",
            "The Rosser-Schoenfeld C=2 entry is a finite-window benchmark, not a global all-range replacement.",
        ],
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
    print("Literature bounds in C/(sqrt(x)log x) scale:")
    for row in rows:
        print(
            f"  x={row.x_prime:>9} "
            f"C_req={row.required_constant:.6f} "
            f"C_actual={row.actual_constant:.6f} "
            f"C=1 margin={row.target_constant_1_margin:+.3f} "
            f"RS margin={row.rosser_schoenfeld_margin:+.3f} "
            f"Dusart C={row.dusart_theorem_6_10_constant:.3f} "
            f"0.2/log^3 C={row.optimistic_0_2_log3_constant:.3f} "
            f"Axler C={row.axler_2018_constant if row.axler_2018_constant is not None else float('nan'):.3f} "
            f"RH C={row.rh_reciprocal_constant_3_over_8pi:.3f}"
        )


if __name__ == "__main__":
    main()
