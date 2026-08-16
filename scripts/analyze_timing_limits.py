"""TESS-only ephemeris and tidal-decay limit for WASP-19 b.

The committed sectors are reduced homogeneously.  Each transit is timed with
the sector-level shape fixed, while midpoint, local baseline, and local slope
remain free.  Time-averaging beta factors and an empirical per-sector timing
jitter are included before comparing linear and quadratic ephemerides.

This is a sensitivity analysis, not a full apsidal-precession fit and not a
replacement for the longer, heterogeneous timing baselines in the literature.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

import analyze_multisector as multi
import analyze_transit as base


TIMINGS_FILE = base.FIG_DIR / "individual_transit_timings.csv"
STATS_FILE = base.FIG_DIR / "timing_limit_statistics.csv"
FIGURE_FILE = base.FIG_DIR / "wasp19b_timing_limits.png"

EARTHS_PER_SUN = 332_946.0
RSUN_PER_AU = 0.00465047
MS_PER_DAY_PER_YEAR = 86_400_000.0 * 365.25


def fit_event(
    time: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    epoch: int,
    sector_shape: np.ndarray,
    sector: int,
) -> dict[str, float | int | bool] | None:
    """Measure one transit midpoint with the sector-level shape held fixed."""
    predicted = base.EPOCH_BJD + epoch * base.PERIOD_DAYS
    offset = time - predicted
    duration = base.DURATION_HOURS / 24.0
    chosen = np.abs(offset) <= 1.8 * duration
    x, y, sigma = offset[chosen], flux[chosen], error[chosen]
    if len(x) < 60:
        return None
    if min(np.sum(x < -0.7 * duration), np.sum(x > 0.7 * duration)) < 10:
        return None
    if np.sum(np.abs(x) < 0.45 * duration) < 10:
        return None

    radius_ratio = float(sector_shape[1])
    impact = float(sector_shape[2])

    def model(parameters: np.ndarray) -> np.ndarray:
        timing, baseline, slope = parameters
        full = np.asarray([timing, radius_ratio, impact, baseline, slope])
        return base.transit_profile(x, full)

    initial = np.asarray([float(sector_shape[0]), float(np.median(y)), 0.0])
    timing_bound = 0.55 * duration
    fit = least_squares(
        lambda pars: (y - model(pars)) / sigma,
        initial,
        bounds=([-timing_bound, 0.94, -0.2], [timing_bound, 1.06, 0.2]),
        x_scale="jac",
        max_nfev=500,
    )
    fitted = model(fit.x)
    null = base.weighted_linear_null(x, y, sigma)
    chi2_transit = float(np.sum(np.square((y - fitted) / sigma)))
    chi2_null = float(np.sum(np.square((y - null) / sigma)))
    dof = len(y) - 3
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
    _, _, _, beta = multi.noise_curve(y - fitted)
    scale = np.sqrt(max(chi2_transit / dof, 1.0)) * beta
    timing_error = float(np.sqrt(max(covariance[0, 0], 0.0)) * scale)
    delta_bic = chi2_null + 2 * np.log(len(y)) - (chi2_transit + 3 * np.log(len(y)))
    return {
        "sector": sector,
        "epoch": epoch,
        "predicted_bjd": predicted,
        "measured_bjd": predicted + float(fit.x[0]),
        "oc_seconds": float(fit.x[0]) * 86_400.0,
        "formal_timing_error_seconds": timing_error * 86_400.0,
        "timing_error_seconds": timing_error * 86_400.0,
        "sector_jitter_seconds": 0.0,
        "n_points": len(y),
        "beta": beta,
        "reduced_chi_square": chi2_transit / dof,
        "delta_bic": delta_bic,
        "supported": bool(delta_bic >= 10 and timing_error > 0),
    }


def weighted_ephemeris(
    epoch: np.ndarray,
    oc_days: np.ndarray,
    error_days: np.ndarray,
    degree: int,
) -> dict[str, object]:
    """Fit a conditioned linear or quadratic O-C ephemeris."""
    center = float(np.round(np.average(epoch, weights=1.0 / error_days**2)))
    x = epoch - center
    columns = [np.ones_like(x), x]
    if degree == 2:
        columns.append(x**2)
    design = np.column_stack(columns)
    weighted_design = design / error_days[:, None]
    coefficients, _, _, _ = np.linalg.lstsq(
        weighted_design, oc_days / error_days, rcond=None
    )
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    model = design @ coefficients
    chi_square = float(np.sum(np.square((oc_days - model) / error_days)))
    return {
        "center_epoch": center,
        "x": x,
        "coefficients": coefficients,
        "covariance": covariance,
        "model": model,
        "chi_square": chi_square,
        "dof": len(epoch) - len(coefficients),
        "bic": chi_square + len(coefficients) * np.log(len(epoch)),
    }


def robust_sigma(values: np.ndarray) -> float:
    median = float(np.median(values))
    return float(1.4826 * np.median(np.abs(values - median)))


def add_sector_jitter(events: list[dict[str, float | int | bool]]) -> dict[int, float]:
    """Inflate errors by the within-sector scatter about a first linear fit."""
    supported = [event for event in events if event["supported"]]
    epoch = np.asarray([event["epoch"] for event in supported], dtype=float)
    oc = np.asarray([event["oc_seconds"] for event in supported], dtype=float) / 86_400.0
    error = np.asarray(
        [event["formal_timing_error_seconds"] for event in supported], dtype=float
    ) / 86_400.0
    first = weighted_ephemeris(epoch, oc, error, degree=1)
    residual_seconds = (oc - np.asarray(first["model"])) * 86_400.0
    jitters: dict[int, float] = {}
    for sector in sorted({int(event["sector"]) for event in supported}):
        selected = np.asarray([int(event["sector"]) == sector for event in supported])
        scatter = robust_sigma(residual_seconds[selected]) if selected.sum() >= 3 else 0.0
        median_formal = float(np.median([
            float(event["formal_timing_error_seconds"])
            for event in supported if int(event["sector"]) == sector
        ]))
        # A small floor prevents a few formally precise events from dominating.
        jitters[sector] = max(scatter, 0.5 * median_formal)
    for event in events:
        jitter = jitters.get(int(event["sector"]), 0.0)
        event["sector_jitter_seconds"] = jitter
        event["timing_error_seconds"] = float(np.hypot(
            float(event["formal_timing_error_seconds"]), jitter
        ))
    return jitters


def tidal_quality_factor(period_dot_days_per_day: float) -> float:
    """Equilibrium-tide Q'_* using the convention stated in the report."""
    with (base.DATA_DIR / "system_parameters.csv").open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle))
    planet_mass_earth = float(row["pl_bmasse"])
    star_mass_sun = float(row["st_mass"])
    star_radius_sun = float(row["st_rad"])
    semimajor_au = float(row["pl_orbsmax"])
    mass_ratio = planet_mass_earth / (star_mass_sun * EARTHS_PER_SUN)
    radius_ratio = star_radius_sun * RSUN_PER_AU / semimajor_au
    numerator = 27.0 * np.pi / 2.0 * mass_ratio * radius_ratio**5
    return float(numerator / abs(period_dot_days_per_day))


def main() -> dict[str, object]:
    base.FIG_DIR.mkdir(exist_ok=True)
    events: list[dict[str, float | int | bool]] = []
    for path in sorted(base.DATA_DIR.glob("tess*_lc.fits"), key=multi.sector_number):
        sector = multi.sector_number(path)
        time, flux, error, _ = base.load_light_curve(path)
        sector_fit = base.compare_models(time, flux, error)
        shape = np.asarray(sector_fit["parameters"], dtype=float)
        first_epoch = int(np.ceil((time.min() - base.EPOCH_BJD) / base.PERIOD_DAYS))
        last_epoch = int(np.floor((time.max() - base.EPOCH_BJD) / base.PERIOD_DAYS))
        for epoch in range(first_epoch, last_epoch + 1):
            event = fit_event(time, flux, error, epoch, shape, sector)
            if event is not None:
                events.append(event)

    supported = [event for event in events if event["supported"]]
    if len(supported) < 12:
        raise RuntimeError("Too few supported transits for the timing comparison")
    jitters = add_sector_jitter(events)

    epoch = np.asarray([event["epoch"] for event in supported], dtype=float)
    oc_days = np.asarray([event["oc_seconds"] for event in supported], dtype=float) / 86_400.0
    error_days = np.asarray(
        [event["timing_error_seconds"] for event in supported], dtype=float
    ) / 86_400.0
    linear = weighted_ephemeris(epoch, oc_days, error_days, degree=1)
    quadratic = weighted_ephemeris(epoch, oc_days, error_days, degree=2)
    coefficient = float(quadratic["coefficients"][2])
    coefficient_error = float(np.sqrt(quadratic["covariance"][2, 2]))
    period_dot = 2.0 * coefficient / base.PERIOD_DAYS
    period_dot_error = 2.0 * coefficient_error / base.PERIOD_DAYS
    period_dot_ms_year = period_dot * MS_PER_DAY_PER_YEAR
    period_dot_error_ms_year = period_dot_error * MS_PER_DAY_PER_YEAR
    negative_95_limit_ms_year = min(period_dot_ms_year - 1.96 * period_dot_error_ms_year, -1e-12)
    negative_95_limit = negative_95_limit_ms_year / MS_PER_DAY_PER_YEAR
    q_lower_95 = tidal_quality_factor(negative_95_limit)

    fields = [
        "sector", "epoch", "predicted_bjd", "measured_bjd", "oc_seconds",
        "formal_timing_error_seconds", "sector_jitter_seconds", "timing_error_seconds",
        "n_points", "beta", "reduced_chi_square", "delta_bic", "supported",
    ]
    with TIMINGS_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(events)

    rows = [
        ("events_fitted", len(events), "count"),
        ("events_supported", len(supported), "count; individual Delta BIC >= 10"),
        ("sectors", len(jitters), "count"),
        ("timing_baseline_days", float(max(event["measured_bjd"] for event in supported) - min(event["measured_bjd"] for event in supported)), "days"),
        ("linear_chi_square", linear["chi_square"], ""),
        ("linear_dof", linear["dof"], ""),
        ("linear_bic", linear["bic"], ""),
        ("quadratic_chi_square", quadratic["chi_square"], ""),
        ("quadratic_dof", quadratic["dof"], ""),
        ("quadratic_bic", quadratic["bic"], ""),
        ("delta_bic_linear_minus_quadratic", linear["bic"] - quadratic["bic"], ""),
        ("period_dot_ms_per_year", period_dot_ms_year, "ms per year; conditional TESS-only estimate"),
        ("period_dot_error_ms_per_year", period_dot_error_ms_year, "ms per year; formal after beta and sector jitter"),
        ("negative_period_dot_95_limit_ms_per_year", negative_95_limit_ms_year, "ms per year; Gaussian conditional bound"),
        ("modified_stellar_tidal_quality_factor_95_lower", q_lower_95, "Q-prime; equilibrium-tide convention"),
    ]
    for sector, jitter in jitters.items():
        rows.append((f"sector_{sector}_timing_jitter_seconds", jitter, "seconds"))
    with STATS_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        for name, value, unit in rows:
            writer.writerow([name, f"{value:.12g}" if isinstance(value, float) else value, unit])

    x = np.asarray(linear["x"])
    grid = np.linspace(x.min() - 20, x.max() + 20, 800)
    linear_curve = linear["coefficients"][0] + linear["coefficients"][1] * grid
    quadratic_curve = (
        quadratic["coefficients"][0]
        + quadratic["coefficients"][1] * grid
        + quadratic["coefficients"][2] * grid**2
    )
    residual_seconds = (oc_days - np.asarray(linear["model"])) * 86_400.0
    fig, (ax, residual_ax) = plt.subplots(2, 1, figsize=(10, 8.2), constrained_layout=True)
    colors = {9: "#2563eb", 62: "#059669", 63: "#d97706"}
    for sector in sorted(jitters):
        selected = np.asarray([int(event["sector"]) == sector for event in supported])
        ax.errorbar(
            x[selected], oc_days[selected] * 86_400.0,
            yerr=error_days[selected] * 86_400.0,
            fmt="o", ms=4, capsize=2, alpha=0.72,
            color=colors.get(sector, "#475569"), label=f"TESS Sector {sector}",
        )
        residual_ax.scatter(
            x[selected], residual_seconds[selected], s=24, alpha=0.72,
            color=colors.get(sector, "#475569"), label=f"S{sector}",
        )
    ax.plot(grid, linear_curve * 86_400.0, "--", color="#475569", lw=1.8, label="linear ephemeris")
    ax.plot(grid, quadratic_curve * 86_400.0, color="#9d174d", lw=2.1, label="quadratic ephemeris")
    ax.set(
        xlabel=f"Transit epoch relative to E = {linear['center_epoch']:.0f}",
        ylabel="Observed - archive prediction [s]",
        title=f"WASP-19 b: {len(supported)} independently timed TESS transits",
    )
    ax.text(
        0.02, 0.04,
        f"Pdot = {period_dot_ms_year:+.2f} +/- {period_dot_error_ms_year:.2f} ms yr^-1\n"
        f"Delta BIC(linear - quadratic) = {linear['bic'] - quadratic['bic']:.1f}\n"
        f"Conditional 95% negative bound = {negative_95_limit_ms_year:.2f} ms yr^-1",
        transform=ax.transAxes, fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cbd5e1"},
    )
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    residual_ax.axhline(0, color="#475569", ls="--", lw=1.3)
    residual_ax.set(
        xlabel=f"Transit epoch relative to E = {linear['center_epoch']:.0f}",
        ylabel="Residual from linear ephemeris [s]",
        title="Residual structure after the preferred linear timing model",
    )
    residual_ax.grid(alpha=0.2)
    fig.savefig(FIGURE_FILE, dpi=190)
    plt.close(fig)

    return {
        "events": events,
        "supported": supported,
        "jitters": jitters,
        "linear": linear,
        "quadratic": quadratic,
        "period_dot_ms_per_year": period_dot_ms_year,
        "period_dot_error_ms_per_year": period_dot_error_ms_year,
        "negative_95_limit_ms_per_year": negative_95_limit_ms_year,
        "q_lower_95": q_lower_95,
    }


if __name__ == "__main__":
    result = main()
    print(
        f"WASP-19 b: {len(result['supported'])} supported transits; "
        f"Pdot={result['period_dot_ms_per_year']:+.3f} +/- "
        f"{result['period_dot_error_ms_per_year']:.3f} ms/yr; "
        f"Delta BIC={result['linear']['bic'] - result['quadratic']['bic']:.2f}"
    )
