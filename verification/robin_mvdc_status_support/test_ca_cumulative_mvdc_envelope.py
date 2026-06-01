from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_CA_CSV = ROOT / "mvdc_ca_beta_target_high_alpha_ca_rows.csv"
DEFAULT_BLOCK_CSV = ROOT / "ca_successive_mvdc_blocks_rows.csv"
DEFAULT_JSON = ROOT / "ca_cumulative_mvdc_envelope_report.json"
DEFAULT_CSV = ROOT / "ca_cumulative_mvdc_envelope_rows.csv"


@dataclass
class CACumulativeEnvelopeRow:
    y_prime: int
    x_prime: int
    beta_error_e_x: float
    ca_reserve_x: float
    previous_upper_e_y: float
    moment_1: float
    m1_step_threshold: float
    m1_step_margin: float
    cumulative_upper_e_x: float
    upper_minus_actual_e: float
    reserve_margin_with_upper: float
    reserve_margin_with_actual_e: float


def read_ca_rows(path: Path) -> list[tuple[int, float, float]]:
    rows: list[tuple[int, float, float]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows.append((int(raw["last_prime"]), float(raw["beta_error_e"]), float(raw["total_a_plus_b"])))
    return rows


def read_m1_blocks(path: Path) -> dict[tuple[int, int], float]:
    blocks: dict[tuple[int, int], float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            blocks[(int(raw["y_prime"]), int(raw["x_prime"]))] = float(raw["moment_1"])
    return blocks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Propagate E(x) upper bounds along sampled CA support using Q<=M1.")
    parser.add_argument("--ca-csv", type=Path, default=DEFAULT_CA_CSV)
    parser.add_argument("--block-csv", type=Path, default=DEFAULT_BLOCK_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ca_rows = read_ca_rows(args.ca_csv)
    blocks = read_m1_blocks(args.block_csv)

    upper = ca_rows[0][1]
    output_rows: list[CACumulativeEnvelopeRow] = []
    for (y, _e_y, _r_y), (x, e_x, r_x) in zip(ca_rows, ca_rows[1:]):
        m1 = blocks[(y, x)]
        previous_upper = upper
        upper = previous_upper + m1
        threshold = r_x - previous_upper
        output_rows.append(
            CACumulativeEnvelopeRow(
                y_prime=y,
                x_prime=x,
                beta_error_e_x=e_x,
                ca_reserve_x=r_x,
                previous_upper_e_y=previous_upper,
                moment_1=m1,
                m1_step_threshold=threshold,
                m1_step_margin=threshold - m1,
                cumulative_upper_e_x=upper,
                upper_minus_actual_e=upper - e_x,
                reserve_margin_with_upper=r_x - upper,
                reserve_margin_with_actual_e=r_x - e_x,
            )
        )

    payload = {
        "notes": [
            "This is the cumulative CA-support version of the safe MVDC envelope.",
            "The only analytic input used per step is Q(Y,x)<=M1(Y,x).",
            "Starting from a certified base E(x0), define U_i=U_{i-1}+M1(x_{i-1},x_i). Then E(x_i)<=U_i.",
            "The reported margin checks whether U_i is still below the CA ledger reserve A_CA+B_log_CA.",
        ],
        "base": {
            "prime": ca_rows[0][0],
            "exact_beta_error_e": ca_rows[0][1],
            "ca_reserve": ca_rows[0][2],
            "base_margin": ca_rows[0][2] - ca_rows[0][1],
        },
        "rows": [asdict(row) for row in output_rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(output_rows[0]).keys()))
        writer.writeheader()
        for row in output_rows:
            writer.writerow(asdict(row))

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Cumulative CA MVDC envelope:")
    for row in output_rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"U={row.cumulative_upper_e_x:+.6e} "
            f"E={row.beta_error_e_x:+.6e} "
            f"M1_margin={row.m1_step_margin:+.3e} "
            f"U-E={row.upper_minus_actual_e:+.3e} "
            f"reserve-U={row.reserve_margin_with_upper:+.6e}"
        )


if __name__ == "__main__":
    main()
