from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from test_mvdc_signed_chebyshev_center import DEFAULT_CSV as DEFAULT_CENTER_CSV
from test_mvdc_signed_chebyshev_ledger import DEFAULT_CSV as DEFAULT_LEDGER_CSV
from test_w_mvdc_certified_envelope import primes_up_to


ROOT = Path(__file__).resolve().parent
DEFAULT_JSON = ROOT / "mvdc_prime_cell_moments_report.json"
DEFAULT_CSV = ROOT / "mvdc_prime_cell_moments_rows.csv"
DEFAULT_DETAIL_CSV = ROOT / "mvdc_prime_cell_moments_detail_rows.csv"


GAUSS_8_NODES = (
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363,
)

GAUSS_8_WEIGHTS = (
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763,
)


@dataclass
class PrimeCellMomentRow:
    y_prime: int
    x_prime: int
    prime_count: int
    required_k_for_ledger: float
    center_k_from_cells: float
    center_k_from_center_csv: float
    k_margin_from_cells: float
    k_margin_from_ledger_csv: float
    total_cell_surplus: float
    normalized_m1_about_required_k: float
    weighted_m2_about_required_k: float
    weighted_m3_about_required_k: float
    weighted_m4_about_required_k: float
    weighted_variance_about_center: float
    weighted_std_about_center: float
    weighted_skewness_about_center: float
    weighted_kurtosis_about_center: float
    min_z: float
    min_z_left_prime: int
    max_z: float
    max_z_left_prime: int
    weighted_q05_z: float
    weighted_q50_z: float
    weighted_q95_z: float
    negative_lambda_mass_about_required_k: float
    positive_lambda_mass_about_required_k: float
    negative_surplus_abs: float
    positive_surplus: float
    min_cell_surplus: float
    min_cell_surplus_left_prime: int
    max_cell_surplus: float
    max_cell_surplus_left_prime: int
    min_partial_surplus: float
    min_partial_surplus_left_prime: int
    final_partial_surplus: float
    min_deficit_d: float
    min_deficit_left_prime: int
    max_deficit_d: float
    max_deficit_left_prime: int
    min_deficit_threshold: float
    min_deficit_threshold_left_prime: int
    max_deficit_threshold: float
    max_deficit_threshold_left_prime: int
    min_deficit_margin: float
    min_deficit_margin_left_prime: int
    abel_initial_component: float
    abel_gap_component: float
    b_weight_sum: float
    required_k_c_weight_sum: float
    abel_surplus: float
    abel_identity_error: float
    c_sum_from_cells: float
    sqrt_weight_integral_from_center_csv: float
    c_sum_relative_error: float
    numerator_sum_from_cells: float
    numerator_sum_from_center_csv: float
    numerator_relative_error: float
    first_moment_identity_error: float


@dataclass
class PrimeCellDetailRow:
    y_prime: int
    x_prime: int
    left_prime: int
    right_prime: int
    theta_left: float
    deficit_d: float
    a_weight: float
    b_weight: float
    c_weight: float
    lambda_weight: float
    z_value: float
    required_k_for_ledger: float
    z_minus_required_k: float
    deficit_threshold: float
    deficit_margin: float
    cell_surplus: float
    partial_surplus: float


def read_center_rows(path: Path) -> dict[tuple[int, int], dict[str, float]]:
    rows: dict[tuple[int, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            y = int(raw["y_prime"])
            x = int(raw["x_prime"])
            rows[(y, x)] = {
                "prime_count": float(raw["prime_count"]),
                "sqrt_weight_integral": float(raw["sqrt_weight_integral"]),
                "mvdc_center_k": float(raw["mvdc_center_k"]),
                "s2_block": float(raw["s2_block"]),
            }
    return rows


def read_ledger_rows(path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows.append(
                {
                    "y_prime": int(raw["y_prime"]),
                    "x_prime": int(raw["x_prime"]),
                    "required_k_for_ledger": float(raw["required_k_for_ledger"]),
                    "k_margin": float(raw["k_margin"]),
                }
            )
    return rows


def prime_indices(primes: list[int]) -> dict[int, int]:
    return {prime: index for index, prime in enumerate(primes)}


def theta_prefixes(primes: Iterable[int]) -> list[float]:
    theta = 0.0
    values: list[float] = []
    for prime in primes:
        theta += math.log(prime)
        values.append(theta)
    return values


def g_weight(t: float) -> float:
    return 1.0 / (t * math.log(t))


def h_weight(t: float) -> float:
    return math.log(math.log(t)) - 1.0 / math.log(t)


def sqrt_weight_density_log_u(u_value: float) -> float:
    return (u_value + 1.0) * math.exp(-0.5 * u_value) / (u_value * u_value)


def c_weight_gauss(left_prime: int, right_prime: int) -> float:
    left = math.log(left_prime)
    right = math.log(right_prime)
    midpoint = 0.5 * (left + right)
    half_width = 0.5 * (right - left)
    total = 0.0
    for node, weight in zip(GAUSS_8_NODES, GAUSS_8_WEIGHTS):
        total += weight * sqrt_weight_density_log_u(midpoint + half_width * node)
    return half_width * total


def weighted_quantiles(values: list[tuple[float, float]], quantiles: tuple[float, ...]) -> list[float]:
    if not values:
        return [math.nan for _ in quantiles]
    ordered = sorted(values, key=lambda item: item[0])
    total_weight = sum(weight for _value, weight in ordered)
    result: list[float] = []
    cumulative = 0.0
    index = 0
    for quantile in quantiles:
        threshold = quantile * total_weight
        while index < len(ordered) - 1 and cumulative + ordered[index][1] < threshold:
            cumulative += ordered[index][1]
            index += 1
        result.append(ordered[index][0])
    return result


def write_summary_csv(path: Path, rows: list[PrimeCellMomentRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def open_detail_writer(path: Path, enabled: bool) -> tuple[csv.DictWriter | None, object | None]:
    if not enabled:
        return None, None
    handle = path.open("w", newline="", encoding="utf-8")
    writer = csv.DictWriter(handle, fieldnames=list(PrimeCellDetailRow.__dataclass_fields__.keys()))
    writer.writeheader()
    return writer, handle


def evaluate_block(
    primes: list[int],
    theta_values: list[float],
    y_index: int,
    x_index: int,
    y: int,
    x: int,
    required_k: float,
    center_row: dict[str, float],
    ledger_margin: float,
    detail_writer: csv.DictWriter | None,
) -> PrimeCellMomentRow:
    c_sum = 0.0
    numerator_sum = 0.0
    partial_surplus = 0.0
    b_weight_sum = 0.0
    abel_gap_component = 0.0
    positive_surplus = 0.0
    negative_surplus_abs = 0.0
    negative_c_weight = 0.0
    positive_c_weight = 0.0

    moment2_required_raw = 0.0
    moment3_required_raw = 0.0
    moment4_required_raw = 0.0
    variance_center_raw = 0.0
    skew_center_raw = 0.0
    kurt_center_raw = 0.0

    min_z = math.inf
    max_z = -math.inf
    min_z_left = 0
    max_z_left = 0
    min_cell_surplus = math.inf
    max_cell_surplus = -math.inf
    min_cell_surplus_left = 0
    max_cell_surplus_left = 0
    min_partial_surplus = math.inf
    min_partial_surplus_left = 0
    min_deficit = math.inf
    max_deficit = -math.inf
    min_deficit_left = 0
    max_deficit_left = 0
    min_deficit_threshold = math.inf
    max_deficit_threshold = -math.inf
    min_deficit_threshold_left = 0
    max_deficit_threshold_left = 0
    min_deficit_margin = math.inf
    min_deficit_margin_left = 0

    center_k_reference = center_row["mvdc_center_k"]
    sqrt_weight_reference = center_row["sqrt_weight_integral"]
    d0 = float(y) - theta_values[y_index]
    g_x = g_weight(x)
    abel_initial_component = d0 * (g_weight(y) - g_x)
    weighted_values: list[tuple[float, float]] = []

    for index in range(y_index, x_index):
        left_prime = primes[index]
        right_prime = primes[index + 1]
        theta_left = theta_values[index]
        deficit_d = float(left_prime) - theta_left

        a_weight = g_weight(left_prime) - g_weight(right_prime)
        b_weight = h_weight(right_prime) - h_weight(left_prime) - float(left_prime) * a_weight
        c_weight = c_weight_gauss(left_prime, right_prime)

        numerator = deficit_d * a_weight + b_weight
        z_value = numerator / c_weight
        z_minus_required = z_value - required_k
        z_minus_center = z_value - center_k_reference
        deficit_threshold = (required_k * c_weight - b_weight) / a_weight
        deficit_margin = deficit_d - deficit_threshold
        cell_surplus = numerator - required_k * c_weight
        partial_surplus += cell_surplus

        c_sum += c_weight
        numerator_sum += numerator
        b_weight_sum += b_weight
        abel_gap_component += (right_prime - left_prime - math.log(right_prime)) * (g_weight(right_prime) - g_x)
        weighted_values.append((z_value, c_weight))
        moment2_required_raw += c_weight * (z_minus_required**2)
        moment3_required_raw += c_weight * (z_minus_required**3)
        moment4_required_raw += c_weight * (z_minus_required**4)
        variance_center_raw += c_weight * (z_minus_center**2)
        skew_center_raw += c_weight * (z_minus_center**3)
        kurt_center_raw += c_weight * (z_minus_center**4)

        if z_value < required_k:
            negative_c_weight += c_weight
        else:
            positive_c_weight += c_weight
        if cell_surplus >= 0.0:
            positive_surplus += cell_surplus
        else:
            negative_surplus_abs += -cell_surplus

        if z_value < min_z:
            min_z = z_value
            min_z_left = left_prime
        if z_value > max_z:
            max_z = z_value
            max_z_left = left_prime
        if cell_surplus < min_cell_surplus:
            min_cell_surplus = cell_surplus
            min_cell_surplus_left = left_prime
        if cell_surplus > max_cell_surplus:
            max_cell_surplus = cell_surplus
            max_cell_surplus_left = left_prime
        if partial_surplus < min_partial_surplus:
            min_partial_surplus = partial_surplus
            min_partial_surplus_left = left_prime
        if deficit_d < min_deficit:
            min_deficit = deficit_d
            min_deficit_left = left_prime
        if deficit_d > max_deficit:
            max_deficit = deficit_d
            max_deficit_left = left_prime
        if deficit_threshold < min_deficit_threshold:
            min_deficit_threshold = deficit_threshold
            min_deficit_threshold_left = left_prime
        if deficit_threshold > max_deficit_threshold:
            max_deficit_threshold = deficit_threshold
            max_deficit_threshold_left = left_prime
        if deficit_margin < min_deficit_margin:
            min_deficit_margin = deficit_margin
            min_deficit_margin_left = left_prime

        if detail_writer is not None:
            detail_writer.writerow(
                asdict(
                    PrimeCellDetailRow(
                    y_prime=y,
                    x_prime=x,
                    left_prime=left_prime,
                    right_prime=right_prime,
                    theta_left=theta_left,
                    deficit_d=deficit_d,
                    a_weight=a_weight,
                    b_weight=b_weight,
                    c_weight=c_weight,
                    lambda_weight=c_weight / sqrt_weight_reference,
                    z_value=z_value,
                    required_k_for_ledger=required_k,
                    z_minus_required_k=z_minus_required,
                    deficit_threshold=deficit_threshold,
                    deficit_margin=deficit_margin,
                    cell_surplus=cell_surplus,
                    partial_surplus=partial_surplus,
                )
                )
            )

    q05, q50, q95 = weighted_quantiles(weighted_values, (0.05, 0.50, 0.95))
    center_k_cells = numerator_sum / c_sum
    normalized_m1 = partial_surplus / c_sum
    variance_center = variance_center_raw / c_sum
    std_center = math.sqrt(max(0.0, variance_center))
    if std_center > 0.0:
        skew_center = (skew_center_raw / c_sum) / (std_center**3)
        kurt_center = (kurt_center_raw / c_sum) / (std_center**4)
    else:
        skew_center = 0.0
        kurt_center = 0.0

    numerator_reference = center_k_reference * sqrt_weight_reference
    required_k_c_weight_sum = required_k * c_sum
    abel_surplus = abel_initial_component + abel_gap_component + b_weight_sum - required_k_c_weight_sum
    abel_identity_error = abel_surplus - partial_surplus
    c_relative_error = (c_sum - sqrt_weight_reference) / sqrt_weight_reference
    numerator_relative_error = (numerator_sum - numerator_reference) / numerator_reference
    first_moment_identity_error = normalized_m1 - (center_k_cells - required_k)

    return PrimeCellMomentRow(
        y_prime=y,
        x_prime=x,
        prime_count=x_index - y_index,
        required_k_for_ledger=required_k,
        center_k_from_cells=center_k_cells,
        center_k_from_center_csv=center_k_reference,
        k_margin_from_cells=center_k_cells - required_k,
        k_margin_from_ledger_csv=ledger_margin,
        total_cell_surplus=partial_surplus,
        normalized_m1_about_required_k=normalized_m1,
        weighted_m2_about_required_k=moment2_required_raw / c_sum,
        weighted_m3_about_required_k=moment3_required_raw / c_sum,
        weighted_m4_about_required_k=moment4_required_raw / c_sum,
        weighted_variance_about_center=variance_center,
        weighted_std_about_center=std_center,
        weighted_skewness_about_center=skew_center,
        weighted_kurtosis_about_center=kurt_center,
        min_z=min_z,
        min_z_left_prime=min_z_left,
        max_z=max_z,
        max_z_left_prime=max_z_left,
        weighted_q05_z=q05,
        weighted_q50_z=q50,
        weighted_q95_z=q95,
        negative_lambda_mass_about_required_k=negative_c_weight / c_sum,
        positive_lambda_mass_about_required_k=positive_c_weight / c_sum,
        negative_surplus_abs=negative_surplus_abs,
        positive_surplus=positive_surplus,
        min_cell_surplus=min_cell_surplus,
        min_cell_surplus_left_prime=min_cell_surplus_left,
        max_cell_surplus=max_cell_surplus,
        max_cell_surplus_left_prime=max_cell_surplus_left,
        min_partial_surplus=min_partial_surplus,
        min_partial_surplus_left_prime=min_partial_surplus_left,
        final_partial_surplus=partial_surplus,
        min_deficit_d=min_deficit,
        min_deficit_left_prime=min_deficit_left,
        max_deficit_d=max_deficit,
        max_deficit_left_prime=max_deficit_left,
        min_deficit_threshold=min_deficit_threshold,
        min_deficit_threshold_left_prime=min_deficit_threshold_left,
        max_deficit_threshold=max_deficit_threshold,
        max_deficit_threshold_left_prime=max_deficit_threshold_left,
        min_deficit_margin=min_deficit_margin,
        min_deficit_margin_left_prime=min_deficit_margin_left,
        abel_initial_component=abel_initial_component,
        abel_gap_component=abel_gap_component,
        b_weight_sum=b_weight_sum,
        required_k_c_weight_sum=required_k_c_weight_sum,
        abel_surplus=abel_surplus,
        abel_identity_error=abel_identity_error,
        c_sum_from_cells=c_sum,
        sqrt_weight_integral_from_center_csv=sqrt_weight_reference,
        c_sum_relative_error=c_relative_error,
        numerator_sum_from_cells=numerator_sum,
        numerator_sum_from_center_csv=numerator_reference,
        numerator_relative_error=numerator_relative_error,
        first_moment_identity_error=first_moment_identity_error,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prime-cell MVDC moment audit for the signed Chebyshev centre K(Y,x). "
            "By default it writes compact block summaries; use --write-detail for one row per prime cell."
        )
    )
    parser.add_argument("--center-csv", type=Path, default=DEFAULT_CENTER_CSV)
    parser.add_argument("--ledger-csv", type=Path, default=DEFAULT_LEDGER_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--detail-csv", type=Path, default=DEFAULT_DETAIL_CSV)
    parser.add_argument("--write-detail", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    center_rows = read_center_rows(args.center_csv)
    ledger_rows = read_ledger_rows(args.ledger_csv)
    if not ledger_rows:
        raise ValueError("No ledger rows found.")

    max_prime = max(int(row["x_prime"]) for row in ledger_rows)
    primes = primes_up_to(max_prime)
    indices = prime_indices(primes)
    theta_values = theta_prefixes(primes)
    detail_writer, detail_handle = open_detail_writer(args.detail_csv, args.write_detail)

    try:
        rows: list[PrimeCellMomentRow] = []
        for ledger_row in ledger_rows:
            y = int(ledger_row["y_prime"])
            x = int(ledger_row["x_prime"])
            key = (y, x)
            if key not in center_rows:
                raise ValueError(f"Missing centre row for block {y}->{x}")
            if y not in indices or x not in indices:
                raise ValueError(f"Block endpoint is not in prime table: {y}->{x}")
            rows.append(
                evaluate_block(
                    primes=primes,
                    theta_values=theta_values,
                    y_index=indices[y],
                    x_index=indices[x],
                    y=y,
                    x=x,
                    required_k=float(ledger_row["required_k_for_ledger"]),
                    center_row=center_rows[key],
                    ledger_margin=float(ledger_row["k_margin"]),
                    detail_writer=detail_writer,
                )
            )
    finally:
        if detail_handle is not None:
            detail_handle.close()

    payload = {
        "notes": [
            "This is the prime-cell MVDC audit requested by Appendix No. 11.",
            "For each cell (q_j,q_{j+1}) it uses A_j, B_j, C_j, D_j=q_j-theta(q_j), Z_j=(D_j A_j+B_j)/C_j.",
            "The first moment around the ledger threshold is K_i-K_i^req and equals total_cell_surplus / sum C_j.",
            "C_j is evaluated by 8-point Gauss-Legendre quadrature in the logarithmic variable; the report checks the sum against the block integral from the centre CSV.",
            "Use --write-detail to emit one CSV row per prime cell; the default output is intentionally compact.",
        ],
        "inputs": {
            "center_csv": str(args.center_csv),
            "ledger_csv": str(args.ledger_csv),
            "write_detail": args.write_detail,
            "detail_csv": str(args.detail_csv) if args.write_detail else None,
        },
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    write_summary_csv(args.csv, rows)

    print(f"primes up to {max_prime}: {len(primes)}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    if args.write_detail:
        print(f"wrote {args.detail_csv}")
    print()
    print("Prime-cell MVDC moment audit:")
    for row in rows:
        print(
            f"{row.y_prime:>9}->{row.x_prime:<9} "
            f"Kcell={row.center_k_from_cells:.9f} "
            f"Kreq={row.required_k_for_ledger:.9f} "
            f"M1={row.normalized_m1_about_required_k:+.9f} "
            f"neg_lambda={row.negative_lambda_mass_about_required_k:.4f} "
            f"minDmargin={row.min_deficit_margin:.3f}@{row.min_deficit_margin_left_prime} "
            f"q05/q50/q95={row.weighted_q05_z:.3f}/{row.weighted_q50_z:.3f}/{row.weighted_q95_z:.3f} "
            f"min_partial={row.min_partial_surplus:+.3e} "
            f"Cerr={row.c_sum_relative_error:+.2e}"
        )


if __name__ == "__main__":
    main()
