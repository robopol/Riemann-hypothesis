#!/usr/bin/env python3
"""Plot the signed prime-power dynamics and the triangular no-go mechanism.

Chart contract
--------------
Question:
    How do the one-sided spectral state, the signed triangular discrepancy,
    and the cone reserves behave at prime-power events through a finite range?
Takeaway:
    The sampled |E_h|/Q envelope narrows strongly, and the finite conditional
    reserve stays positive after a small base, but a growing exponential mode
    can still have a vanishing relative triangular remainder.
Surface:
    Reproducible static Matplotlib figure (PNG and SVG) plus sampled CSV data.
Grain:
    Logarithmically sampled prime-power events, augmented by all early events;
    exact decade envelopes are read from the archived full binary64 scan.

The chart is exploratory.  It is not interval arithmetic and does not prove
an asymptotic envelope, the universal Turan inequality, or RH.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from exploratory_eventwise_reserve_scan import (
    EULER_GAMMA,
    LOG_TWO_PI,
    WEIL_C,
    trivial_f_array,
)
from exploratory_signed_triangular_scan import (
    build_event_support,
    lambda_key,
    phi_at_x,
    required_support_limit,
)


DEFAULT_REPORT = Path(__file__).with_name(
    "exploratory_signed_triangular_scan_1e8_report.json"
)
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if (
    SCRIPT_DIRECTORY.name == "finite_robin_prime_power"
    and SCRIPT_DIRECTORY.parent.name == "verification"
):
    DEFAULT_OUTPUT_DIR = SCRIPT_DIRECTORY.parents[1] / "papers" / "figures"
else:
    DEFAULT_OUTPUT_DIR = SCRIPT_DIRECTORY.parent / "output"

BLUE = "#2F5DA8"
ORANGE = "#D97824"
GOLD = "#A67C00"
INK = "#20252B"
NEUTRAL = "#66717E"
LIGHT_NEUTRAL = "#AAB2BC"
GRID = "#DDE1E6"


def sample_event_indices(
    event_q: np.ndarray,
    event_exponent: np.ndarray,
    center_count: int,
    limit: int,
    sample_count: int,
    early_limit: int = 10_000,
) -> np.ndarray:
    """Return log-spaced indices plus all early and higher-power events."""

    targets = np.geomspace(2.0, float(limit), sample_count)
    indices = np.searchsorted(event_q[:center_count], targets, side="left")
    indices = np.clip(indices, 0, center_count - 1)
    previous = np.maximum(indices - 1, 0)
    use_previous = np.abs(np.log(event_q[previous]) - np.log(targets)) < np.abs(
        np.log(event_q[indices]) - np.log(targets)
    )
    indices = np.where(use_previous, previous, indices)
    early_count = int(
        np.searchsorted(event_q[:center_count], early_limit, side="right")
    )
    early = np.arange(early_count, dtype=np.int64)
    higher_powers = np.flatnonzero(event_exponent[:center_count] >= 2)
    return np.unique(
        np.concatenate((early, higher_powers, indices.astype(np.int64)))
    )


def evaluate_sample(
    limit: int,
    lambda_value: float,
    sample_count: int,
    segment_span: int,
) -> dict[str, np.ndarray]:
    """Evaluate the signed triangular state on sampled prime-power events."""

    support_limit = required_support_limit(limit, lambda_value)
    event_q, event_prime, event_exponent, event_lambda = build_event_support(
        support_limit, segment_span
    )
    cumulative_psi = np.cumsum(event_lambda, dtype=np.float64)
    cumulative_s_lambda = np.cumsum(
        event_lambda / event_q.astype(np.float64), dtype=np.float64
    )
    center_count = int(np.searchsorted(event_q, limit, side="right"))
    indices = sample_event_indices(
        event_q, event_exponent, center_count, limit, sample_count
    )

    q = event_q[indices]
    prime = event_prime[indices]
    exponent = event_exponent[indices]
    lambda_q = event_lambda[indices]
    x = q.astype(np.float64)
    root_x = np.sqrt(x)
    log_x = np.log(x)
    y = 1.0 / x

    psi_right = cumulative_psi[indices]
    s_right = cumulative_s_lambda[indices]
    f_value = trivial_f_array(y)
    f_prime = np.log1p(-(y * y))
    u_right = (
        s_right
        - log_x
        + EULER_GAMMA
        + 0.5 * y * f_prime
        - 0.5 * f_value
    )
    v_right = psi_right - x + LOG_TWO_PI + 0.5 * f_prime
    spectral_b = root_x * u_right - v_right / root_x
    phi = WEIL_C - spectral_b
    p_right = -0.5 * (root_x * u_right + v_right / root_x)
    jump = lambda_q / root_x
    p_bar = p_right + 0.5 * jump

    b_factor = 2.0 * (1.0 - 1.0 / (x * (x * x - 1.0)))
    radius_squared = phi * (
        root_x * b_factor - 0.5 * WEIL_C + 0.25 * phi
    )
    if np.any(radius_squared <= 0.0):
        raise RuntimeError("The sampled cone radius is not positive")
    radius = np.sqrt(radius_squared)

    h = lambda_value * log_x / root_x
    phi_minus = phi_at_x(
        x * np.exp(-h), event_q, cumulative_psi, cumulative_s_lambda
    )
    phi_plus = phi_at_x(
        x * np.exp(h), event_q, cumulative_psi, cumulative_s_lambda
    )
    phi_2h = phi_at_x(
        np.exp(2.0 * h), event_q, cumulative_psi, cumulative_s_lambda
    )
    difference = phi_plus - phi_minus
    d_h = difference / (2.0 * h)
    e_h = p_bar - d_h
    turan_bound = np.sqrt(phi * phi_2h) / h
    turan_defect = 4.0 * phi * phi_2h - difference * difference
    conditional_reserve = radius - turan_bound + e_h - 0.5 * jump
    actual_post_reserve = radius + p_right

    return {
        "q": q,
        "prime": prime,
        "exponent": exponent,
        "t": log_x,
        "B": spectral_b,
        "P": phi,
        "Q": radius,
        "E_over_Q": e_h / radius,
        "S_over_Q": turan_bound / radius,
        "M_over_Q": conditional_reserve / radius,
        "actual_reserve_over_Q": actual_post_reserve / radius,
        "normalized_U": turan_defect / (4.0 * phi * phi_2h),
    }


def load_decade_envelope(
    report_path: Path, lambda_value: float
) -> dict[str, np.ndarray | float]:
    """Load exact decade envelopes, the finite fit, and defect minimum."""

    report = json.loads(report_path.read_text(encoding="utf-8"))
    result = report["lambda_results"][lambda_key(lambda_value)]
    centers = []
    maxima = []
    central_minima = []
    reserve_minima = []
    for row in result["decade_diagnostics"]:
        left = max(2, int(row["q_left"]))
        right = int(row["q_right"])
        centers.append(math.sqrt(left * right))
        maxima.append(row["extrema"]["maximum_abs_E_h_over_R"]["value"])
        central_minima.append(
            row["extrema"]["minimum_relative_R_minus_turan_bound"]["value"]
        )
        reserve_minima.append(
            row["extrema"]["minimum_M_T_over_R"]["value"]
        )
    minimum_normalized_u = result["global_extrema"][
        "minimum_normalized_U"
    ]["value"]
    fit = next(
        row
        for row in result["effective_scaling_fits"]
        if row["metric"] == "maximum_abs_E_h_over_R"
        and int(row["cutoff"]) == 10_000
    )
    return {
        "q": np.asarray(centers, dtype=np.float64),
        "maximum_abs_E_over_Q": np.asarray(maxima, dtype=np.float64),
        "minimum_central_allowance_over_Q": np.asarray(
            central_minima, dtype=np.float64
        ),
        "minimum_M_over_Q": np.asarray(reserve_minima, dtype=np.float64),
        "minimum_normalized_U": float(minimum_normalized_u),
        "fit_slope": float(fit["slope_log_abs_E_over_R_on_log_q"]),
        "fit_intercept": float(fit["intercept"]),
        "fit_R_squared": float(fit["R_squared"]),
    }


def stable_relative_triangular_remainder(z: np.ndarray) -> np.ndarray:
    """Return |sinh(z)/z-1| without cancellation near zero."""

    absolute_z = np.abs(z)
    series = z * z / 6.0 + z**4 / 120.0 + z**6 / 5040.0
    direct = np.sinh(z) / z - 1.0
    return np.where(absolute_z < 0.05, np.abs(series), np.abs(direct))


def write_sample_csv(path: Path, data: dict[str, np.ndarray]) -> None:
    """Write the plotted sample with explicit column names."""

    columns = list(data)
    row_count = len(data[columns[0]])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for index in range(row_count):
            writer.writerow([data[column][index] for column in columns])


def style_axis(axis: plt.Axes) -> None:
    """Apply the shared research-chart axis style."""

    axis.set_facecolor("white")
    axis.grid(True, which="major", color=GRID, linewidth=0.7, alpha=0.8)
    axis.grid(True, which="minor", color=GRID, linewidth=0.35, alpha=0.35)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(LIGHT_NEUTRAL)
    axis.spines["bottom"].set_color(LIGHT_NEUTRAL)
    axis.tick_params(colors=INK, labelsize=9)


def plot_figure(
    data: dict[str, np.ndarray],
    diagnostics: dict[str, np.ndarray | float],
    lambda_value: float,
    output_png: Path,
    output_svg: Path,
) -> None:
    """Render the four-panel signed-dynamics figure."""

    q = data["q"].astype(np.float64)
    t = data["t"]
    decade_q = np.asarray(diagnostics["q"])
    decade_envelope = np.asarray(diagnostics["maximum_abs_E_over_Q"])
    minimum_normalized_u = float(diagnostics["minimum_normalized_U"])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "text.color": INK,
            "legend.frameon": False,
            "svg.hashsalt": "riemann-signed-triangular-2026-08-05",
        }
    )
    figure, axes = plt.subplots(2, 2, figsize=(15.5, 10.5), dpi=180)
    figure.patch.set_facecolor("white")

    for axis in axes.flat:
        style_axis(axis)
        axis.set_xscale("log")
        axis.set_xlim(2, float(np.max(q)))

    top_left, top_right, bottom_left, bottom_right = axes.flat

    top_left.plot(
        q,
        data["B"],
        color=ORANGE,
        linewidth=0.75,
        alpha=0.85,
        label=r"$\mathfrak{B}(\log q)$",
    )
    top_left.plot(
        q,
        data["P"],
        color=BLUE,
        linewidth=0.75,
        alpha=0.85,
        label=r"$\mathfrak{P}(\log q)=C-\mathfrak{B}(\log q)$",
    )
    top_left.axhline(0.0, color=INK, linewidth=0.8)
    top_left.set_title("Spectral residual and one-sided state", loc="left")
    top_left.set_ylabel("state value")
    top_left.legend(loc="upper right", fontsize=9)

    positive_e = data["E_over_Q"] >= 0.0
    higher_power = data["exponent"] >= 2
    top_right.scatter(
        q[positive_e],
        np.abs(data["E_over_Q"][positive_e]),
        s=2.0,
        color=BLUE,
        alpha=0.20,
        linewidths=0,
        label=r"sampled $|E_h|/Q$, $E_h>0$",
    )
    top_right.scatter(
        q[~positive_e],
        np.abs(data["E_over_Q"][~positive_e]),
        s=2.0,
        color=ORANGE,
        alpha=0.18,
        linewidths=0,
        label=r"sampled $|E_h|/Q$, $E_h<0$",
    )
    top_right.scatter(
        q[higher_power],
        np.abs(data["E_over_Q"][higher_power]),
        s=15.0,
        facecolors="none",
        edgecolors=GOLD,
        alpha=0.55,
        linewidths=0.55,
        label=r"higher prime powers $p^a$, $a\geq2$",
    )
    top_right.plot(
        decade_q,
        decade_envelope,
        color=ORANGE,
        linewidth=1.8,
        marker="o",
        markersize=3.5,
        label=r"full-scan decade max $|E_h|/Q$",
    )
    fit_q = np.geomspace(10_000.0, float(np.max(q)), 300)
    fitted_envelope = np.exp(float(diagnostics["fit_intercept"])) * fit_q ** float(
        diagnostics["fit_slope"]
    )
    top_right.plot(
        fit_q,
        fitted_envelope,
        color=INK,
        linewidth=1.2,
        linestyle=":",
        label=(
            rf"descriptive fit $q^{{{float(diagnostics['fit_slope']):.3f}}}$, "
            + rf"$R^2={float(diagnostics['fit_R_squared']):.3f}$"
        ),
    )
    top_right.set_yscale("log")
    top_right.set_ylim(1e-8, 1.0)
    top_right.set_title("Magnitude of the signed triangular discrepancy", loc="left")
    top_right.set_ylabel(r"$|E_h|/Q$ at $\lambda=1/4$ (log scale)")
    top_right.legend(loc="lower left", fontsize=7.8, markerscale=2.0)

    bottom_left.scatter(
        q,
        data["actual_reserve_over_Q"],
        s=2.0,
        color=BLUE,
        alpha=0.22,
        linewidths=0,
        label=r"actual $(Q+p^+)/Q$",
    )
    bottom_left.scatter(
        q,
        data["M_over_Q"],
        s=2.0,
        color=ORANGE,
        alpha=0.20,
        linewidths=0,
        label=r"conditional $M_\lambda/Q$",
    )
    bottom_left.scatter(
        q,
        1.0 - data["S_over_Q"],
        s=2.0,
        color=NEUTRAL,
        alpha=0.18,
        linewidths=0,
        label=r"central allowance $(Q-S_h)/Q$",
    )
    bottom_left.plot(
        decade_q,
        np.asarray(diagnostics["minimum_M_over_Q"]),
        color=ORANGE,
        marker="o",
        markersize=3.5,
        linewidth=1.5,
        label=r"full-scan decade min $M_\lambda/Q$",
    )
    bottom_left.plot(
        decade_q,
        np.asarray(diagnostics["minimum_central_allowance_over_Q"]),
        color=NEUTRAL,
        marker="s",
        markersize=3.2,
        linewidth=1.2,
        linestyle="--",
        label=r"full-scan decade min $(Q-S_h)/Q$",
    )
    bottom_left.axhline(0.0, color=INK, linewidth=0.9)
    bottom_left.set_ylim(-1.0, 1.15)
    bottom_left.set_title("Normalized lower-cone reserves", loc="left")
    bottom_left.set_ylabel("fraction of cone radius")
    bottom_left.legend(loc="lower right", fontsize=7.8, markerscale=2.0)
    bottom_left.text(
        0.02,
        0.05,
        (
            r"Finite scan: $M_{1/4}>0$ from $q=23$; "
            + rf"minimum normalized $U={minimum_normalized_u:.3f}$"
        ),
        transform=bottom_left.transAxes,
        fontsize=8.5,
        color=INK,
        bbox={"facecolor": "white", "edgecolor": GRID, "alpha": 0.9},
    )

    mu = 0.25
    h_mode = lambda_value * t * np.exp(-0.5 * t)
    relative_remainder = stable_relative_triangular_remainder(mu * h_mode)
    growing_mode = np.exp(mu * (t - t[0]))
    bottom_right.plot(
        q,
        growing_mode,
        color=BLUE,
        linewidth=1.5,
        label=rf"growing mode $e^{{{mu:g}t}}$ (normalized)",
    )
    bottom_right.plot(
        q,
        relative_remainder,
        color=ORANGE,
        linewidth=1.5,
        linestyle="--",
        label=r"relative triangular remainder",
    )
    bottom_right.set_yscale("log")
    bottom_right.set_title(
        "Why remainder decay alone is insufficient", loc="left"
    )
    bottom_right.set_ylabel("dimensionless magnitude (log scale)")
    bottom_right.legend(loc="center left", fontsize=8.5)

    for axis in axes[1, :]:
        axis.set_xlabel(r"prime-power event $q$ (log scale)")

    figure.suptitle(
        r"Signed prime-power dynamics through $10^8$"
        + rf"  ($h={lambda_value:g}\,\log q/\sqrt{{q}}$)",
        fontsize=18,
        fontweight="bold",
        x=0.06,
        ha="left",
        y=0.985,
    )
    figure.text(
        0.06,
        0.947,
        (
            "Binary64 exploratory data. Log-spaced event sample; orange "
            "envelope uses every event in the archived full scan."
        ),
        fontsize=10.5,
        color=NEUTRAL,
        ha="left",
    )
    figure.text(
        0.06,
        0.015,
        (
            "Finite narrowing is descriptive, not an asymptotic theorem. "
            "The conditional reserve uses the Turan secant scale; its "
            "universal prime-side proof remains open."
        ),
        fontsize=9.5,
        color=NEUTRAL,
        ha="left",
    )
    figure.subplots_adjust(
        left=0.07, right=0.985, top=0.91, bottom=0.075, hspace=0.28, wspace=0.20
    )
    figure.savefig(output_png, dpi=220, facecolor="white")
    figure.savefig(
        output_svg,
        facecolor="white",
        metadata={"Date": None, "Creator": "plot_signed_triangular_dynamics.py"},
    )
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=100_000_000)
    parser.add_argument("--lambda-value", type=float, default=0.25)
    parser.add_argument("--samples", type=int, default=30_000)
    parser.add_argument("--segment-span", type=int, default=10_000_000)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    """Evaluate the sample, export it, and render the figure."""

    args = parse_args()
    if args.limit < 2 or args.samples < 100:
        raise SystemExit("--limit must be >=2 and --samples must be >=100")
    if not (0.0 < args.lambda_value < math.sqrt(2.0)):
        raise SystemExit("--lambda-value must satisfy 0<lambda<sqrt(2)")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = evaluate_sample(
        args.limit,
        args.lambda_value,
        args.samples,
        args.segment_span,
    )
    diagnostics = load_decade_envelope(args.report, args.lambda_value)

    lambda_slug = str(args.lambda_value).replace(".", "p")
    stem = f"signed_triangular_dynamics_lambda_{lambda_slug}"
    output_png = args.output_dir / f"{stem}.png"
    output_svg = args.output_dir / f"{stem}.svg"
    output_csv = args.output_dir / f"{stem}_sample.csv"
    write_sample_csv(output_csv, data)
    plot_figure(
        data,
        diagnostics,
        args.lambda_value,
        output_png,
        output_svg,
    )
    print(f"wrote {output_png}")
    print(f"wrote {output_svg}")
    print(f"wrote {output_csv} ({len(data['q']):,} sampled events)")


if __name__ == "__main__":
    main()
