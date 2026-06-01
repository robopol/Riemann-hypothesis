from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_PROFILE_CSV = ROOT / "mvdc_beta_requirement_profiles.csv"
DEFAULT_CA_CSV = ROOT / "ca_like_ledger_rows.csv"
DEFAULT_OPTIMAL_CSV = ROOT / "optimal_ledger_envelope_rows.csv"
DEFAULT_JSON = ROOT / "beta_env_gap_absorption_report.json"
DEFAULT_CSV = ROOT / "beta_env_gap_absorption_rows.csv"

EULER_GAMMA = 0.577215664901532860606512090082402431
E_GAMMA = math.exp(EULER_GAMMA)


@dataclass
class GapAbsorptionRow:
    source: str
    family: str
    label: str
    x: int
    log_n: float
    bridge_valid: bool
    beta_error_e: float
    eta_label: str
    eta_value: float
    profile_deficit_a: float
    bridge_b: float
    total_reserve: float
    log_gate_margin: float
    beta_plus: float
    envelope_plus: float
    bridge_target: float
    beta_surplus: float
    envelope_gap: float
    additive_margin: float
    identity_residual: float
    multiplicative_margin: float


def bridge_reserve(log_n: float, x: int) -> float:
    if log_n <= x:
        return 0.0
    return math.log(math.log(log_n) / math.log(x))


def beta_plus(x: int, eta: float) -> float:
    return E_GAMMA * math.log(x) * math.exp(eta)


def rosser_eta(x: int, c_m: float) -> float:
    log_x = math.log(x)
    return math.log1p(c_m / (log_x * log_x))


def build_gap_row(
    *,
    source: str,
    family: str,
    label: str,
    x: int,
    log_n: float,
    beta_error_e: float,
    eta_label: str,
    eta_value: float,
    profile_deficit_a: float,
    bridge_b: float,
) -> GapAbsorptionRow:
    total_reserve = profile_deficit_a + bridge_b
    log_gate_margin = total_reserve - eta_value
    beta_upper = beta_plus(x, eta_value)
    envelope_upper = beta_upper * math.exp(-profile_deficit_a)
    bridge_target = E_GAMMA * math.log(x) * math.exp(bridge_b)
    beta_surplus = beta_upper - bridge_target
    envelope_gap = beta_upper - envelope_upper
    additive_margin = envelope_gap - beta_surplus
    identity_rhs = E_GAMMA * math.log(x) * math.exp(eta_value - profile_deficit_a)
    identity_rhs *= math.expm1(log_gate_margin)
    identity_residual = additive_margin - identity_rhs
    multiplicative_margin = bridge_target / envelope_upper - 1.0
    return GapAbsorptionRow(
        source=source,
        family=family,
        label=label,
        x=x,
        log_n=log_n,
        bridge_valid=log_n > x,
        beta_error_e=beta_error_e,
        eta_label=eta_label,
        eta_value=eta_value,
        profile_deficit_a=profile_deficit_a,
        bridge_b=bridge_b,
        total_reserve=total_reserve,
        log_gate_margin=log_gate_margin,
        beta_plus=beta_upper,
        envelope_plus=envelope_upper,
        bridge_target=bridge_target,
        beta_surplus=beta_surplus,
        envelope_gap=envelope_gap,
        additive_margin=additive_margin,
        identity_residual=identity_residual,
        multiplicative_margin=multiplicative_margin,
    )


def eta_values(x: int, beta_error_e: float) -> list[tuple[str, float]]:
    return [
        ("actual_E_plus", max(beta_error_e, 0.0)),
        ("rosser_half", rosser_eta(x, 0.5)),
        ("rosser_one", rosser_eta(x, 1.0)),
    ]


def read_profile_rows(path: Path) -> list[GapAbsorptionRow]:
    rows: list[GapAbsorptionRow] = []
    if not path.exists():
        return rows
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            family = raw["family"]
            if family not in {"floor_log_lcm", "ca_like"}:
                continue
            x = int(raw["last_prime"])
            log_n = float(raw["log_n"])
            beta_error_e = float(raw["beta_error_E"])
            deficit = float(raw["true_deficit_A"])
            bridge = bridge_reserve(log_n, x)
            for eta_label, eta_value in eta_values(x, beta_error_e):
                rows.append(
                    build_gap_row(
                        source=path.name,
                        family=family,
                        label=raw["label"],
                        x=x,
                        log_n=log_n,
                        beta_error_e=beta_error_e,
                        eta_label=eta_label,
                        eta_value=eta_value,
                        profile_deficit_a=deficit,
                        bridge_b=bridge,
                    )
                )
    return rows


def read_ca_rows(path: Path) -> list[GapAbsorptionRow]:
    rows: list[GapAbsorptionRow] = []
    if not path.exists():
        return rows
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            x = int(raw["last_prime"])
            log_n = x + float(raw["log_n_minus_last_prime"])
            beta_error_e = float(raw["beta_error_e"])
            deficit = float(raw["true_deficit_a"])
            bridge = float(raw["bridge_b"])
            for eta_label, eta_value in eta_values(x, beta_error_e):
                rows.append(
                    build_gap_row(
                        source=path.name,
                        family="ca_like",
                        label=f"eps={raw['epsilon']}",
                        x=x,
                        log_n=log_n,
                        beta_error_e=beta_error_e,
                        eta_label=eta_label,
                        eta_value=eta_value,
                        profile_deficit_a=deficit,
                        bridge_b=bridge,
                    )
                )
    return rows


def read_optimal_rows(path: Path) -> list[GapAbsorptionRow]:
    rows: list[GapAbsorptionRow] = []
    if not path.exists():
        return rows
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            x = int(raw["x"])
            beta_error_e = float(raw["beta_error_e"])
            variants = [
                (
                    "continuous_optimal",
                    "water_filling",
                    float(raw["optimal_log_n"]),
                    float(raw["optimal_deficit_a"]),
                    float(raw["optimal_bridge_b"]),
                ),
                (
                    "floor_log",
                    "floor_log_lcm",
                    float("nan"),
                    float(raw["floor_log_deficit_a"]),
                    float(raw["floor_log_bridge_b"]),
                ),
            ]
            for family, label, log_n, deficit, bridge in variants:
                if math.isnan(log_n):
                    # Reconstruct only the bridge target from bridge_b.
                    log_n = math.exp(math.log(x) * math.exp(bridge))
                for eta_label, eta_value in eta_values(x, beta_error_e):
                    rows.append(
                        build_gap_row(
                            source=path.name,
                            family=family,
                            label=label,
                            x=x,
                            log_n=log_n,
                            beta_error_e=beta_error_e,
                            eta_label=eta_label,
                            eta_value=eta_value,
                            profile_deficit_a=deficit,
                            bridge_b=bridge,
                        )
                    )
    return rows


def write_csv(path: Path, rows: list[GapAbsorptionRow]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit beta^+ - B_env^+ absorption for certified full-profile envelopes."
    )
    parser.add_argument("--profile-csv", type=Path, default=DEFAULT_PROFILE_CSV)
    parser.add_argument("--ca-csv", type=Path, default=DEFAULT_CA_CSV)
    parser.add_argument("--optimal-csv", type=Path, default=DEFAULT_OPTIMAL_CSV)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = [
        *read_profile_rows(args.profile_csv),
        *read_ca_rows(args.ca_csv),
        *read_optimal_rows(args.optimal_csv),
    ]
    write_csv(args.csv, rows)
    payload = {
        "notes": [
            "log_gate_margin = A_env + B_log_env - eta. Nonnegative closes the logarithmic gate.",
            "additive_margin = (beta^+ - B_env^+) - (beta^+ - e^gamma log L_env). It is equivalent to the log gate.",
            "identity_residual checks additive_margin = e^gamma log(x) exp(eta-A) expm1(log_gate_margin).",
            "actual_E_plus is the ideal sharp beta envelope eta=max(E,0); Rosser rows show why the old surplus is too coarse.",
            "continuous_optimal rows are relaxed lower-envelope diagnostics; negative margins indicate that this model is too adversarial or needs stronger structural restrictions.",
        ],
        "max_abs_identity_residual": max(
            (abs(row.identity_residual) for row in rows),
            default=0.0,
        ),
        "rows": [asdict(row) for row in rows],
    }
    args.json.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    print()
    print("Beta-envelope gap absorption, selected actual_E_plus rows:")
    selected = [
        row
        for row in rows
        if row.eta_label == "actual_E_plus"
        and row.family in {"ca_like", "floor_log_lcm", "continuous_optimal", "floor_log"}
    ]
    for row in selected[:80]:
        print(
            f"{row.source:34s} {row.family:18s} {row.label:14s} "
            f"x={row.x:<9} A+B-eta={row.log_gate_margin:+.3e} "
            f"gap-surplus={row.additive_margin:+.3e} bridge={row.bridge_valid}"
        )


if __name__ == "__main__":
    main()
