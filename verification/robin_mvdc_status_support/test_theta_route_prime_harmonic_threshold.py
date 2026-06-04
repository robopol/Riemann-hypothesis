from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_TREND_JSON = ROOT / "c2_shortfall_trend_report.json"
DEFAULT_TREND_CSV = ROOT / "c2_shortfall_trend_rows.csv"
DEFAULT_JSON = ROOT / "theta_route_prime_harmonic_threshold_report.json"
DEFAULT_CSV = ROOT / "theta_route_prime_harmonic_threshold_rows.csv"


@dataclass
class ThetaRoutePrimeHarmonicThresholdRow:
    y_prime: int
    x_prime: int
    block_prime_count: int
    route_mu: float
    scale_sqrt_x_log_x: float
    actual_constant: float
    c_req_actual_reserve: float
    ca_reserve_actual: float
    ca_deficit_actual: float
    ca_bridge_log: float
    theta_phi: float
    theta_phi_over_deficit: float
    c_req_theta: float
    theta_margin: float
    theta_drop_from_actual_reserve: float
    passes_theta_threshold: bool


def dusart_theta_delta(t: float) -> float:
    if t <= 1.0:
        raise ValueError("Theta bounds are used only for t > 1.")

    log_t = math.log(t)
    candidates = [1.0]
    if t > 2.0:
        candidates.append(1.2323 / log_t)
        candidates.append(3.965 / (log_t * log_t))
    if t > 908_994.0:
        candidates.append(0.001 / log_t)
    if t > 3_594_641.0:
        candidates.append(0.2 / (log_t * log_t))
    if t > 122_568_683.0:
        candidates.append(0.05 / (log_t * log_t))
    if t > 7_713_133_853.0:
        candidates.append(0.01 / (log_t * log_t))
    return min(candidates)


def theta_lower(t: float) -> float:
    return t * (1.0 - dusart_theta_delta(t))


def theta_upper(t: float) -> float:
    return t * (1.0 + dusart_theta_delta(t))


def theta_support_upper(x: int) -> float:
    log_x = math.log(x)
    return x * (1.0 + 1.0 / (25.0 * log_x * log_x))


def layer_endpoint(x: int, layer: int) -> float:
    if layer == 1:
        return float(x)
    if layer < 1:
        raise ValueError("layer must be positive.")

    x_plus = theta_support_upper(x)
    target = (x_plus + 1.0) * math.log(x_plus)
    low = 1.0 + 1e-12
    high = float(x)
    for _ in range(180):
        mid = (low + high) / 2.0
        value = (mid**layer) * math.log(mid)
        if value < target:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def theta_block_layer_lower(a: float, b: float, layer: int, ratio: float) -> float:
    if b <= a:
        return 0.0
    if ratio <= 1.0:
        raise ValueError("--theta-block-ratio must be greater than 1.")

    total = 0.0
    left = a
    guard = 0
    while left < b:
        right = min(b, ratio * left)
        theta_mass_lower = max(theta_lower(right) - theta_upper(left), 0.0)
        total += theta_mass_lower / (math.log(right) * (right**layer))
        if right == b:
            break
        if right <= left:
            raise RuntimeError("Theta block loop stopped making progress.")
        left = right
        guard += 1
        if guard > 100_000:
            raise RuntimeError("Too many theta blocks; check the block ratio.")
    return total


def theta_phi(x: int, layer_order: int, ratio: float) -> float:
    endpoints = {layer: layer_endpoint(x, layer) for layer in range(1, layer_order + 1)}
    return sum(
        theta_block_layer_lower(endpoints[layer], endpoints[layer - 1], layer, ratio)
        for layer in range(2, layer_order + 1)
    )


def ca_bridge_from_point(point: dict[str, object]) -> float:
    x = int(point["last_prime"])
    log_n_minus_x = float(point["log_n_minus_last_prime"])
    log_n = x + log_n_minus_x
    if log_n <= x:
        return 0.0
    return math.log(math.log(log_n) / math.log(x))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Substitute the finite Dusart theta-block CA-deficit envelope into "
            "the Route-II reciprocal-prime harmonic threshold."
        )
    )
    parser.add_argument("--trend-json", type=Path, default=DEFAULT_TREND_JSON)
    parser.add_argument("--trend-csv", type=Path, default=DEFAULT_TREND_CSV)
    parser.add_argument("--layer-order", type=int, default=6)
    parser.add_argument("--theta-block-ratio", type=float, default=1.30)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.layer_order < 2:
        raise ValueError("--layer-order must be at least 2.")

    trend_report = json.loads(args.trend_json.read_text(encoding="utf-8"))
    ca_points = {
        int(point["last_prime"]): point
        for point in trend_report.get("ca_points", [])
    }
    with args.trend_csv.open(newline="", encoding="utf-8") as handle:
        trend_rows = list(csv.DictReader(handle))
    if not trend_rows:
        raise ValueError(f"No rows found in {args.trend_csv}")

    output_rows: list[ThetaRoutePrimeHarmonicThresholdRow] = []
    for raw in trend_rows:
        y = int(raw["y_prime"])
        x = int(raw["x_prime"])
        point = ca_points.get(x)
        if point is None:
            raise KeyError(f"Missing CA point for x={x} in {args.trend_json}")

        block_prime_count = int(raw["block_prime_count"])
        route_mu = math.log(math.log(x) / math.log(y)) / block_prime_count
        scale = math.sqrt(x) * math.log(x)
        actual_constant = float(raw["actual_constant"])
        c_req_actual = float(raw["required_constant"])

        ca_reserve_actual = float(point["reserve_a_plus_b"])
        ca_bridge_log = ca_bridge_from_point(point)
        ca_deficit_actual = ca_reserve_actual - ca_bridge_log
        phi = theta_phi(x, args.layer_order, args.theta_block_ratio)
        c_req_theta = c_req_actual + math.exp(route_mu) * (phi - ca_deficit_actual) * scale
        theta_margin = c_req_theta - actual_constant

        output_rows.append(
            ThetaRoutePrimeHarmonicThresholdRow(
                y_prime=y,
                x_prime=x,
                block_prime_count=block_prime_count,
                route_mu=route_mu,
                scale_sqrt_x_log_x=scale,
                actual_constant=actual_constant,
                c_req_actual_reserve=c_req_actual,
                ca_reserve_actual=ca_reserve_actual,
                ca_deficit_actual=ca_deficit_actual,
                ca_bridge_log=ca_bridge_log,
                theta_phi=phi,
                theta_phi_over_deficit=phi / ca_deficit_actual,
                c_req_theta=c_req_theta,
                theta_margin=theta_margin,
                theta_drop_from_actual_reserve=c_req_actual - c_req_theta,
                passes_theta_threshold=theta_margin >= 0.0,
            )
        )

    payload = {
        "notes": [
            "This is a numerical Route-II substitution audit, not an infinite-range proof.",
            "Input c2_shortfall_trend rows provide C_req from the actual CA reserve.",
            "Input c2_shortfall_trend report provides CA support points with reserve_a_plus_b and log_n_minus_last_prime.",
            "The script splits reserve_a_plus_b into ca_deficit_actual + ca_bridge_log.",
            "theta_phi is the finite prime-free Dusart theta-block lower envelope for ca_deficit_actual.",
            "c_req_theta = c_req_actual_reserve + exp(mu)*(theta_phi-ca_deficit_actual)*sqrt(x)*log(x).",
            "passes_theta_threshold means C_actual <= C_req_theta on the sampled CA-support block.",
        ],
        "args": {
            "trend_json": str(args.trend_json),
            "trend_csv": str(args.trend_csv),
            "layer_order": args.layer_order,
            "theta_block_ratio": args.theta_block_ratio,
        },
        "summary": {
            "rows": len(output_rows),
            "failures": sum(not row.passes_theta_threshold for row in output_rows),
            "min_theta_margin": min(row.theta_margin for row in output_rows),
            "max_theta_margin": max(row.theta_margin for row in output_rows),
            "min_c_req_theta": min(row.c_req_theta for row in output_rows),
            "max_c_req_theta": max(row.c_req_theta for row in output_rows),
            "min_theta_drop": min(row.theta_drop_from_actual_reserve for row in output_rows),
            "max_theta_drop": max(row.theta_drop_from_actual_reserve for row in output_rows),
            "min_theta_phi_over_deficit": min(row.theta_phi_over_deficit for row in output_rows),
            "max_theta_phi_over_deficit": max(row.theta_phi_over_deficit for row in output_rows),
        },
        "rows": [asdict(row) for row in output_rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(output_rows[0]).keys()))
        writer.writeheader()
        for row in output_rows:
            writer.writerow(asdict(row))

    print(f"read {args.trend_json}")
    print(f"read {args.trend_csv}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Theta Route-II prime-harmonic threshold:")
    for row in output_rows:
        print(
            f"  {row.y_prime:>9}->{row.x_prime:<9} "
            f"C_req_act={row.c_req_actual_reserve:.6f} "
            f"C_req_theta={row.c_req_theta:.6f} "
            f"C_actual={row.actual_constant:.6f} "
            f"margin={row.theta_margin:+.6f} "
            f"phi/D={row.theta_phi_over_deficit:.4f} "
            f"pass={row.passes_theta_threshold}"
        )


if __name__ == "__main__":
    main()
