from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_timing_limits as timing


def test_weighted_ephemeris_recovers_quadratic_term():
    epoch = np.arange(-30, 31, dtype=float)
    error = np.full_like(epoch, 2e-6)
    truth = 1e-4 + 3e-7 * epoch - 7e-10 * epoch**2
    fit = timing.weighted_ephemeris(epoch, truth, error, degree=2)
    assert np.allclose(fit["coefficients"], [1e-4, 3e-7, -7e-10], atol=1e-14)


def test_tidal_quality_factor_is_positive():
    assert 1e4 < timing.tidal_quality_factor(-1e-9) < 1e8


def test_real_tess_timing_analysis_writes_auditable_outputs():
    result = timing.main()
    assert len(result["supported"]) >= 30
    assert len(result["jitters"]) == 3
    assert np.isfinite(result["period_dot_ms_per_year"])
    assert result["period_dot_error_ms_per_year"] > 0
    assert result["q_lower_95"] > 0
    assert timing.TIMINGS_FILE.stat().st_size > 1000
    assert timing.STATS_FILE.stat().st_size > 200
    assert timing.FIGURE_FILE.stat().st_size > 10_000
