from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_THETA_CSV = ROOT / "ca_theta_cancellation_rows.csv"
DEFAULT_CENTER_CSV = ROOT / "mvdc_signed_chebyshev_center_rows.csv"
DEFAULT_JSON = ROOT / "mvdc_signed_chebyshev_ledger_report.json"
DEFAULT_CSV = ROOT / "mvdc_signed_chebyshev_ledger_rows.csv"


@dataclass
class SignedChebyshevLedgerRow:
    y_prime: int
    x_prime: int
    previous_upper_w: float
    reserve_x: float
    mvdc_center_k: float
    required_k_for_ledger: float
    k_margin: float
    sqrt_weight_integral: float
    s2_block: float
    corrected_block_from_k: float
    cumulative_upper_w: float
    exact_w_x: float
    upper_minus_exact_w: float
    reserve_margin_with_upper: float
    reserve_margin_with_exact_w: float


def read_theta_rows(path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            prime = int(raw["last_prime"])
            rows[prime] = {key: float(value) for key, value in raw.items() if key != "last_prime"}
    return rows


def read_center_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "y_prime": int(raw["y_prime"]),
                    "x_prime": int(raw["x_prime"]),
                    "mvdc_center_k": float(raw["mvdc_center_k"]),
                    "sqrt_weight_integral": float(raw["sqrt_weight_integral"]),
                    "s2_block": float(raw["s2_block"]),
                }
            )
    return rows


def write_csv(path: Path, rows: list[SignedChebyshevLedgerRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cumulative CA ledger written directly in MVDC signed-Chebyshev centre form.")
    parser.add_argument("--theta-csv", type=Path, default=DEFAULT_THETA_CSV)
    parser.add_argument("--center-csv", type=Path, default=DEFAULT_CENTER_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    theta_rows = read_theta_rows(args.theta_csv)
    center_rows = read_center_rows(args.center_csv)
    if not center_rows:
        raise ValueError("No MVDC signed Chebyshev centre rows found.")

    first_prime = int(center_rows[0]["y_prime"])
    upper = theta_rows[first_prime]["modified_error_w"]
    rows: list[SignedChebyshevLedgerRow] = []

    for raw in center_rows:
        y = int(raw["y_prime"])
        x = int(raw["x_prime"])
        center_k = float(raw["mvdc_center_k"])
        sqrt_weight = float(raw["sqrt_weight_integral"])
        s2_value = float(raw["s2_block"])
        reserve_x = theta_rows[x]["reserve_after_cancellation"]
        exact_w_x = theta_rows[x]["modified_error_w"]

        required_k = (upper + s2_value - reserve_x) / sqrt_weight
        corrected_block = -center_k * sqrt_weight + s2_value
        previous_upper = upper
        upper = previous_upper + corrected_block

        rows.append(
            SignedChebyshevLedgerRow(
                y_prime=y,
                x_prime=x,
                previous_upper_w=previous_upper,
                reserve_x=reserve_x,
                mvdc_center_k=center_k,
                required_k_for_ledger=required_k,
                k_margin=center_k - required_k,
                sqrt_weight_integral=sqrt_weight,
                s2_block=s2_value,
                corrected_block_from_k=corrected_block,
                cumulative_upper_w=upper,
                exact_w_x=exact_w_x,
                upper_minus_exact_w=upper - exact_w_x,
                reserve_margin_with_upper=reserve_x - upper,
                reserve_margin_with_exact_w=reserve_x - exact_w_x,
            )
        )

    payload = {
        "notes": [
            "This is the cumulative CA ledger in the signed Chebyshev MVDC centre form.",
            "The corrected block is -K_i * W_i + S2_i, where K_i is the MVDC centre of t-theta(t) in sqrt(t) units.",
            "The proof target on each step is K_i >= (U_{i-1}+S2_i-R_i)/W_i.",
            "This avoids the artificial fixed C=2 local-step target and uses the actual accumulated CA reserve.",
        ],
        "base": {
            "prime": first_prime,
            "exact_w": theta_rows[first_prime]["modified_error_w"],
            "reserve": theta_rows[first_prime]["reserve_after_cancellation"],
            "base_margin": theta_rows[first_prime]["reserve_after_cancellation"] - theta_rows[first_prime]["modified_error_w"],
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_csv(args.csv, rows)

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("MVDC signed Chebyshev cumulative ledger:")
    for row in rows:
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"K={row.mvdc_center_k:.6f} KreqLedger={row.required_k_for_ledger:.6f} "
            f"Kmargin={row.k_margin:+.6f} "
            f"R-U={row.reserve_margin_with_upper:+.6e} "
            f"U-W={row.upper_minus_exact_w:+.2e}"
        )


if __name__ == "__main__":
    main()
