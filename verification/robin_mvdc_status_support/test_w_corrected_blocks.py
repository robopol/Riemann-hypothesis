from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_THETA_CSV = ROOT / "ca_theta_cancellation_rows.csv"
DEFAULT_BLOCK_CSV = ROOT / "ca_successive_mvdc_blocks_rows.csv"
DEFAULT_JSON = ROOT / "w_corrected_blocks_report.json"
DEFAULT_CSV = ROOT / "w_corrected_blocks_rows.csv"


@dataclass
class WCorrectedBlockRow:
    y_prime: int
    x_prime: int
    w_y: float
    w_x: float
    q_exact: float
    moment_1: float
    delta_rho: float
    w_block_exact: float
    w_block_m1_upper: float
    identity_error: float
    c_step_exact: float
    c_step_m1: float
    scaled_w_y: float
    scaled_w_x: float


def read_theta_rows(path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            x = int(raw["last_prime"])
            rows[x] = {key: float(value) for key, value in raw.items() if key != "last_prime"}
    return rows


def read_block_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "y_prime": int(raw["y_prime"]),
                    "x_prime": int(raw["x_prime"]),
                    "q_exact": float(raw["q_exact"]),
                    "moment_1": float(raw["moment_1"]),
                }
            )
    return rows


def f_scale(x: int) -> float:
    return 1.0 / (math.sqrt(x) * math.log(x))


def scaled_w(x: int, w_value: float) -> float:
    return w_value / f_scale(x)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check corrected W-block descent on sampled CA support.")
    parser.add_argument("--theta-csv", type=Path, default=DEFAULT_THETA_CSV)
    parser.add_argument("--block-csv", type=Path, default=DEFAULT_BLOCK_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    theta_rows = read_theta_rows(args.theta_csv)
    block_rows = read_block_rows(args.block_csv)

    rows: list[WCorrectedBlockRow] = []
    for block in block_rows:
        y = int(block["y_prime"])
        x = int(block["x_prime"])
        if y not in theta_rows or x not in theta_rows:
            continue
        y_row = theta_rows[y]
        x_row = theta_rows[x]
        w_y = y_row["modified_error_w"]
        w_x = x_row["modified_error_w"]
        rho_y = y_row["endpoint_term"]
        rho_x = x_row["endpoint_term"]
        q_exact = float(block["q_exact"])
        moment_1 = float(block["moment_1"])
        delta_rho = rho_x - rho_y
        w_block_exact = q_exact - delta_rho
        w_block_m1_upper = moment_1 - delta_rho
        denominator = f_scale(y) - f_scale(x)
        rows.append(
            WCorrectedBlockRow(
                y_prime=y,
                x_prime=x,
                w_y=w_y,
                w_x=w_x,
                q_exact=q_exact,
                moment_1=moment_1,
                delta_rho=delta_rho,
                w_block_exact=w_block_exact,
                w_block_m1_upper=w_block_m1_upper,
                identity_error=w_block_exact - (w_x - w_y),
                c_step_exact=-w_block_exact / denominator,
                c_step_m1=-w_block_m1_upper / denominator,
                scaled_w_y=scaled_w(y, w_y),
                scaled_w_x=scaled_w(x, w_x),
            )
        )

    payload = {
        "notes": [
            "W=E-(theta(x)-x)/(x log x).",
            "The exact corrected block is W(x)-W(y)=Q(y,x)-(rho(x)-rho(y)).",
            "The safe MVDC upper block replaces Q by M1.",
            "c_step_m1 is the C for which M1-delta_rho <= -C(F(y)-F(x)).",
        ],
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    if rows:
        with args.csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
            writer.writeheader()
            for row in rows:
                writer.writerow(asdict(row))

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Corrected W blocks:")
    for row in rows:
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"Wscale {row.scaled_w_y:.3f}->{row.scaled_w_x:.3f} "
            f"C_exact={row.c_step_exact:.3f} "
            f"C_m1={row.c_step_m1:.3f} "
            f"id_err={row.identity_error:+.2e}"
        )


if __name__ == "__main__":
    main()
