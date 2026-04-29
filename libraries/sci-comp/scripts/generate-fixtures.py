"""Generate JSON fixtures for the TypeScript stats port.

Mirrors the test cases in
  c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/

For each validate_*.py we emit one JSON file:
  src/stats/__tests__/fixtures/<name>.json

Reference values are computed live from scipy / numpy / a self-coded Fisher
oracle. The output is fully self-describing — see README.md.

Run with the SEND-TEST conda env:
  $USERPROFILE/anaconda3/envs/SEND-TEST/python.exe scripts/generate-fixtures.py
"""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from math import comb
from pathlib import Path
from typing import Any

import numpy as np
from scipy import stats as sp_stats


HERE = Path(__file__).resolve().parent
LIB_ROOT = HERE.parent
FIXTURES_DIR = LIB_ROOT / "src" / "stats" / "__tests__" / "fixtures"


# ── Tolerances (mirror validate_helpers.py) ───────────────────────
TOL_PVALUE = 1e-3
TOL_STATISTIC = 1e-2


# ══════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════

def _ndarray_to_list(x: Any) -> Any:
    """Recursive: numpy arrays / scalars → Python lists / floats."""
    if isinstance(x, np.ndarray):
        return [_ndarray_to_list(v) for v in x.tolist()]
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, dict):
        return {k: _ndarray_to_list(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_ndarray_to_list(v) for v in x]
    if isinstance(x, float) and math.isnan(x):
        return None  # serialize NaN as null; TS test re-interprets
    return x


def _meta(source: str, tolerances: dict[str, float] | None = None) -> dict:
    return {
        "source": source,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "scipy_version": sp_stats.__name__ and __import__("scipy").__version__,
        "tolerances": tolerances or {
            "p_value": TOL_PVALUE,
            "statistic": TOL_STATISTIC,
        },
    }


def _write(name: str, data: dict) -> None:
    path = FIXTURES_DIR / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_ndarray_to_list(data), f, indent=2, allow_nan=False, ensure_ascii=False)
    n_cases = len(data.get("cases", []))
    print(f"  wrote {path.name:40s} ({n_cases} cases)")


# ══════════════════════════════════════════════════════════════════
# Reference datasets (shared across multiple fixtures)
# ══════════════════════════════════════════════════════════════════

# ── validate_rest_statistics.py ──
TEMP_2PM = [69, 70, 66, 63, 68, 70, 69, 67, 62, 63, 76, 59, 62, 62, 75, 62, 72, 63]
TEMP_5PM = [68, 62, 67, 68, 69, 67, 61, 59, 62, 61, 69, 66, 62, 62, 61, 70]
ROSETTA_D1 = [27.5, 21.0, 19.0, 23.6, 17.0, 17.9, 16.9, 20.1, 21.9, 22.6, 23.1, 19.6, 19.0, 21.7, 21.4]
ROSETTA_D2 = [27.1, 22.0, 20.8, 23.4, 23.4, 23.5, 25.8, 22.0, 24.8, 20.2, 21.9, 22.1, 22.9, 20.5, 24.4]
ROSETTA_D7 = [30.02, 29.99, 30.11, 29.97, 30.01, 29.99]
ROSETTA_D8 = [29.89, 29.93, 29.72, 29.98, 30.02, 29.98]
ROSETTA_X = [3.0, 4.0, 1.0, 2.1]
ROSETTA_Y = [490.2, 340.0, 433.9]
SHIER_MALES = [19, 22, 16, 29, 24]
SHIER_FEMALES = [20, 11, 17, 12]
STAT_DRUG = [3, 5, 1, 4, 3, 5]
STAT_PLACEBO = [4, 8, 6, 2, 1, 9]
AB94_CAREER = [4, 10, 3, 1, 9, 2, 6, 7, 8, 5]
AB94_PSYCH = [5, 8, 6, 2, 10, 3, 9, 4, 7, 1]
WIKI_IQ = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110]
WIKI_TV = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17]
ABDI_RAW_P = [0.000040, 0.016100, 0.612300]
GARCIA_RAW_P = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074]

# ── validate_hedges_g.py: Fisher's Iris ──
IRIS_SETOSA_Y1 = [
    5.1, 4.9, 4.7, 4.6, 5.0, 5.4, 4.6, 5.0, 4.4, 4.9,
    5.4, 4.8, 4.8, 4.3, 5.8, 5.7, 5.4, 5.1, 5.7, 5.1,
    5.4, 5.1, 4.6, 5.1, 4.8, 5.0, 5.0, 5.2, 5.2, 4.7,
    4.8, 5.4, 5.2, 5.5, 4.9, 5.0, 5.5, 4.9, 4.4, 5.1,
    5.0, 4.5, 4.4, 5.0, 5.1, 4.8, 5.1, 4.6, 5.3, 5.0,
]
IRIS_SETOSA_Y2 = [
    3.5, 3.0, 3.2, 3.1, 3.6, 3.9, 3.4, 3.4, 2.9, 3.1,
    3.7, 3.4, 3.0, 3.0, 4.0, 4.4, 3.9, 3.5, 3.8, 3.8,
    3.4, 3.7, 3.6, 3.3, 3.4, 3.0, 3.4, 3.5, 3.4, 3.2,
    3.1, 3.4, 4.1, 4.2, 3.1, 3.2, 3.5, 3.6, 3.0, 3.4,
    3.5, 2.3, 3.2, 3.5, 3.8, 3.0, 3.8, 3.2, 3.7, 3.3,
]
IRIS_VERSI_Y1 = [
    7.0, 6.4, 6.9, 5.5, 6.5, 5.7, 6.3, 4.9, 6.6, 5.2,
    5.0, 5.9, 6.0, 6.1, 5.6, 6.7, 5.6, 5.8, 6.2, 5.6,
    5.9, 6.1, 6.3, 6.1, 6.4, 6.6, 6.8, 6.7, 6.0, 5.7,
    5.5, 5.5, 5.8, 6.0, 5.4, 6.0, 6.7, 6.3, 5.6, 5.5,
    5.5, 6.1, 5.8, 5.0, 5.6, 5.7, 5.7, 6.2, 5.1, 5.7,
]
IRIS_VERSI_Y2 = [
    3.2, 3.2, 3.1, 2.3, 2.8, 2.8, 3.3, 2.4, 2.9, 2.7,
    2.0, 3.0, 2.2, 2.9, 2.9, 3.1, 3.0, 2.7, 2.2, 2.5,
    3.2, 2.8, 2.5, 2.8, 2.9, 3.0, 2.8, 3.0, 2.9, 2.6,
    2.4, 2.4, 2.7, 2.7, 3.0, 3.4, 3.1, 2.3, 3.0, 2.5,
    2.6, 3.0, 2.6, 2.3, 2.7, 3.0, 2.9, 2.9, 2.5, 2.8,
]
IRIS_VIRG_Y1 = [
    6.3, 5.8, 7.1, 6.3, 6.5, 7.6, 4.9, 7.3, 6.7, 7.2,
    6.5, 6.4, 6.8, 5.7, 5.8, 6.4, 6.5, 7.7, 7.7, 6.0,
    6.9, 5.6, 7.7, 6.3, 6.7, 7.2, 6.2, 6.1, 6.4, 7.2,
    7.4, 7.9, 6.4, 6.3, 6.1, 7.7, 6.3, 6.4, 6.0, 6.9,
    6.7, 6.9, 5.8, 6.8, 6.7, 6.7, 6.3, 6.5, 6.2, 5.9,
]
IRIS_VIRG_Y2 = [
    3.3, 2.7, 3.0, 2.9, 3.0, 3.0, 2.5, 2.9, 2.5, 3.6,
    3.2, 2.7, 3.0, 2.5, 2.8, 3.2, 3.0, 3.8, 2.6, 2.2,
    3.2, 2.8, 2.8, 2.7, 3.3, 3.2, 2.8, 3.0, 2.8, 3.0,
    2.8, 3.8, 2.8, 2.8, 2.6, 3.0, 3.4, 3.1, 3.0, 3.1,
    3.1, 3.1, 2.7, 3.2, 3.3, 3.0, 2.5, 3.0, 3.4, 3.0,
]

# ── validate_trend_test_incidence*.py: published cases ──
TANG_COUNTS = [0, 0, 1, 3]
TANG_TOTALS = [12, 12, 12, 12]
DDT4_COUNTS = [4, 1, 1, 7]
DDT4_TOTALS = [56, 49, 48, 41]
DDT_COUNTS = [4, 1, 1, 7, 17, 12, 15]
DDT_TOTALS = [56, 49, 48, 41, 47, 44, 50]
YOUNG_A = ([0, 0, 10], [50, 50, 50])
YOUNG_B = ([0, 5, 10], [50, 50, 50])
YOUNG_C = ([0, 10, 10], [50, 50, 50])

# ── validate_ancova.py: Montgomery Table 15.10 ──
MONTGOMERY_STR = [36, 41, 39, 42, 49, 40, 48, 39, 45, 44, 35, 37, 42, 34, 32]
MONTGOMERY_DIA = [20, 25, 24, 25, 32, 22, 28, 22, 30, 28, 21, 23, 26, 21, 15]
MONTGOMERY_GRP = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]

ESS_X = [
    8.2, 9.4, 7.7, 8.5, 8.2, 6.0, 9.1, 10.1, 6.8, 7.0, 9.7, 9.9,
    5.7, 5.5, 10.2, 10.3, 6.1, 7.0, 8.7, 8.1, 7.6, 10.1, 9.0, 10.5,
]
ESS_Y = [
    287, 290, 254, 307, 271, 209, 243, 348, 234, 210, 286, 371,
    189, 205, 312, 375, 210, 276, 279, 344, 222, 301, 238, 357,
]
ESS_GRP = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]


# ══════════════════════════════════════════════════════════════════
# Generators — one per fixture file
# ══════════════════════════════════════════════════════════════════

def gen_welch_t() -> None:
    """validate_rest_statistics.py — Test groups 1-3."""
    cases: list[dict] = []

    # Group 1: published references
    pub_cases = [
        ("[McDonald] body temperatures 2pm vs 5pm", TEMP_2PM, TEMP_5PM, {"t": 1.3109, "p": 0.1995, "tol": 1e-3}),
        ("[Rosetta] d1 vs d2", ROSETTA_D1, ROSETTA_D2, {"p": 0.021378001462867, "tol_p": 1e-6}),
        ("[Rosetta] d7 vs d8", ROSETTA_D7, ROSETTA_D8, {"p": 0.090773324285671, "tol_p": 1e-6}),
        ("[Rosetta] x vs y (extreme variance ratio)", ROSETTA_X, ROSETTA_Y, {"t": -9.5595, "p": 0.010751561149785, "tol_t": 0.01, "tol_p": 1e-6}),
    ]
    for name, x, y, ref in pub_cases:
        t, p = sp_stats.ttest_ind(x, y, equal_var=False)
        cases.append({
            "name": name,
            "category": "published_reference",
            "inputs": {"x": x, "y": y},
            "expected": {"statistic": float(t), "p_value": float(p)},
            "published": ref,
        })

    # Group 2: scipy cross-validation (basic, hand-calc, unequal sizes, numpy input)
    cv_cases = [
        ("scipy: basic two-group", [10.2, 11.1, 9.8, 10.5, 10.0], [15.3, 14.8, 16.1, 15.0, 15.5]),
        ("hand-calc: t = 3/√2", [4.0, 6.0], [1.0, 3.0]),
        ("scipy: unequal sizes", [1.0, 2.0, 3.0], [10.0, 11.0, 12.0, 13.0, 14.0, 15.0]),
        ("scipy: numpy-style input", [3.1, 4.1, 5.9, 2.6], [7.0, 8.0, 9.0, 10.0]),
    ]
    for name, x, y in cv_cases:
        t, p = sp_stats.ttest_ind(x, y, equal_var=False)
        cases.append({
            "name": name,
            "category": "scipy_crosscheck",
            "inputs": {"x": x, "y": y},
            "expected": {"statistic": float(t), "p_value": float(p)},
        })

    # Group 3: NaN / edge cases (the spec — what should the impl do)
    edge_cases = [
        {
            "name": "NaN in group1 stripped",
            "category": "nan_handling",
            "inputs": {"x": [1.0, None, 3.0, 2.0], "y": [10.0, 11.0, 12.0]},
            "expected_clean_inputs": {"x": [1.0, 3.0, 2.0], "y": [10.0, 11.0, 12.0]},
        },
        {
            "name": "NaN in group2 stripped",
            "category": "nan_handling",
            "inputs": {"x": [10.0, 11.0, 12.0], "y": [None, 1.0, None, 3.0, 2.0]},
            "expected_clean_inputs": {"x": [10.0, 11.0, 12.0], "y": [1.0, 3.0, 2.0]},
        },
    ]
    for ec in edge_cases:
        ci = ec["expected_clean_inputs"]
        t, p = sp_stats.ttest_ind(ci["x"], ci["y"], equal_var=False)
        ec["expected"] = {"statistic": float(t), "p_value": float(p)}
        cases.append(ec)

    # None-returning edge cases
    for nm, x, y in [
        ("single element group1 → None", [5.0], [1.0, 2.0, 3.0]),
        ("single element group2 → None", [1.0, 2.0, 3.0], [5.0]),
        ("empty group → None", [], [1.0, 2.0]),
        ("all-NaN group → None", [None, None], [1.0, 2.0, 3.0]),
    ]:
        cases.append({
            "name": nm,
            "category": "edge_none",
            "inputs": {"x": x, "y": y},
            "expected": {"statistic": None, "p_value": None},
        })

    # Zero variance → NaN (scipy behavior)
    cases.append({
        "name": "zero variance (constant groups) → NaN",
        "category": "edge_nan",
        "inputs": {"x": [5.0, 5.0, 5.0], "y": [5.0, 5.0, 5.0]},
        "expected": {"statistic": "NaN", "p_value": "NaN"},
    })

    _write("welch-t", {"metadata": _meta("validate_rest_statistics.py — Welch t"), "cases": cases})


def gen_mann_whitney() -> None:
    """validate_rest_statistics.py — Test groups 4-5."""
    cases: list[dict] = []

    # Group 4: published references
    pub = [
        ("[Shier] SPSS: diabetes age", SHIER_MALES, SHIER_FEMALES),
        ("[Statology] tied ranks: drug vs placebo", STAT_DRUG, STAT_PLACEBO),
        ("fully separated [1,2,3] vs [4,5,6]", [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]),
    ]
    for name, x, y in pub:
        u, p = sp_stats.mannwhitneyu(x, y, alternative="two-sided")
        cases.append({
            "name": name,
            "category": "published_reference",
            "inputs": {"x": x, "y": y, "alternative": "two-sided"},
            "expected": {"statistic": float(u), "p_value": float(p)},
        })

    # Group 5: NaN handling and edge cases
    cases.append({
        "name": "NaN removal (pairwise)",
        "category": "nan_handling",
        "inputs": {"x": [1.0, None, 2.0, 3.0], "y": [10.0, 11.0, None]},
        "expected_clean_inputs": {"x": [1.0, 2.0, 3.0], "y": [10.0, 11.0]},
        "expected": (lambda: (lambda u, p: {"statistic": float(u), "p_value": float(p)})(*sp_stats.mannwhitneyu([1.0, 2.0, 3.0], [10.0, 11.0], alternative="two-sided")))(),
    })
    for nm, x, y in [
        ("empty group → None", [], [1.0, 2.0]),
        ("all-NaN group → None", [None], [1.0, 2.0]),
    ]:
        cases.append({
            "name": nm,
            "category": "edge_none",
            "inputs": {"x": x, "y": y, "alternative": "two-sided"},
            "expected": {"statistic": None, "p_value": None},
        })

    _write("mann-whitney", {"metadata": _meta("validate_rest_statistics.py — Mann-Whitney U"), "cases": cases})


def gen_spearman() -> None:
    """validate_rest_statistics.py — Test groups 6-7."""
    cases: list[dict] = []

    # Group 6: published references
    rho_ab94, p_ab94 = sp_stats.spearmanr(AB94_CAREER, AB94_PSYCH)
    cases.append({
        "name": "[AB94] p.466: clinical psychology rankings (n=10)",
        "category": "published_reference",
        "inputs": {"x": AB94_CAREER, "y": AB94_PSYCH},
        "expected": {"rho": float(rho_ab94), "p_value": float(p_ab94)},
        "published": {"rho_exact_rational": "113/165", "rho_decimal": 113.0/165.0, "tol": 1e-10, "sum_d_squared": 52.0, "t_approx": 2.658},
    })

    rho_wiki, p_wiki = sp_stats.spearmanr(WIKI_IQ, WIKI_TV)
    cases.append({
        "name": "[WikiSpear]: IQ vs TV (negative correlation)",
        "category": "published_reference",
        "inputs": {"x": WIKI_IQ, "y": WIKI_TV},
        "expected": {"rho": float(rho_wiki), "p_value": float(p_wiki)},
        "published": {"rho_exact_rational": "-29/165", "rho_decimal": -29.0/165.0, "tol": 1e-10, "p_threshold_max": 0.05},
    })

    # Perfect ±1
    for nm, x, y, exp_rho in [
        ("perfect positive", [1, 2, 3, 4, 5], [10, 20, 30, 40, 50], 1.0),
        ("perfect negative", [1, 2, 3, 4, 5], [50, 40, 30, 20, 10], -1.0),
    ]:
        rho, p = sp_stats.spearmanr(x, y)
        cases.append({
            "name": nm,
            "category": "perfect",
            "inputs": {"x": x, "y": y},
            "expected": {"rho": float(rho), "p_value": float(p) if not math.isnan(p) else "NaN"},
            "expected_rho_exact": exp_rho,
        })

    # Group 7: NaN / tied / edge
    rho_clean, p_clean = sp_stats.spearmanr([1.0, 4.0, 5.0], [10.0, 40.0, 50.0])
    cases.append({
        "name": "NaN pairwise removal",
        "category": "nan_handling",
        "inputs": {"x": [1.0, None, 3.0, 4.0, 5.0], "y": [10.0, 20.0, None, 40.0, 50.0]},
        "expected_clean_inputs": {"x": [1.0, 4.0, 5.0], "y": [10.0, 40.0, 50.0]},
        "expected": {"rho": float(rho_clean), "p_value": float(p_clean) if not math.isnan(p_clean) else "NaN"},
    })

    rho_tied, p_tied = sp_stats.spearmanr([1, 2, 2, 3, 4], [5, 5, 6, 7, 8])
    cases.append({
        "name": "tied ranks",
        "category": "tied",
        "inputs": {"x": [1, 2, 2, 3, 4], "y": [5, 5, 6, 7, 8]},
        "expected": {"rho": float(rho_tied), "p_value": float(p_tied)},
    })

    for nm, x, y in [
        ("fewer than 3 pairs → None", [1, 2], [3, 4]),
        ("all NaN → None", [None], [None]),
    ]:
        cases.append({
            "name": nm,
            "category": "edge_none",
            "inputs": {"x": x, "y": y},
            "expected": {"rho": None, "p_value": None},
        })

    _write("spearman", {"metadata": _meta("validate_rest_statistics.py — Spearman correlation"), "cases": cases})


def gen_severity_trend() -> None:
    """validate_rest_statistics.py — Test group 8 (Spearman wrapper with constant-y guard)."""
    cases: list[dict] = []

    for nm, dl, sev in [
        ("perfect increasing", [0, 10, 50, 100], [0.1, 0.5, 1.2, 2.0]),
        ("perfect decreasing", [0, 10, 50, 100], [2.0, 1.5, 1.0, 0.5]),
        ("arbitrary [0,1,5,10,50] vs [...]", [0, 1, 5, 10, 50], [0.2, 0.3, 0.8, 0.5, 1.5]),
        ("toxicology dose-response (positive plateau)", [0, 50, 100, 200, 500], [0.0, 0.2, 0.8, 1.5, 1.6]),
    ]:
        rho, p = sp_stats.spearmanr(dl, sev)
        cases.append({
            "name": nm,
            "category": "severity",
            "inputs": {"dose_levels": dl, "avg_severities": sev},
            "expected": {"rho": float(rho), "p_value": float(p)},
        })

    cases.append({
        "name": "constant severity → None (correlation undefined)",
        "category": "edge_constant",
        "inputs": {"dose_levels": [0, 10, 50, 100], "avg_severities": [1.0, 1.0, 1.0, 1.0]},
        "expected": {"rho": None, "p_value": None},
    })

    rho_n, p_n = sp_stats.spearmanr([0, 100, 200], [0.1, 1.5, 2.0])
    cases.append({
        "name": "NaN removal",
        "category": "nan_handling",
        "inputs": {"dose_levels": [0, None, 50, 100, 200], "avg_severities": [0.1, 0.5, None, 1.5, 2.0]},
        "expected_clean_inputs": {"dose_levels": [0, 100, 200], "avg_severities": [0.1, 1.5, 2.0]},
        "expected": {"rho": float(rho_n), "p_value": float(p_n)},
    })

    cases.append({
        "name": "fewer than 3 pairs → None",
        "category": "edge_none",
        "inputs": {"dose_levels": [0, 10], "avg_severities": [0.5, 1.0]},
        "expected": {"rho": None, "p_value": None},
    })

    _write("severity-trend", {"metadata": _meta("validate_rest_statistics.py — severity_trend"), "cases": cases})


def gen_welch_pairwise() -> None:
    """validate_rest_statistics.py — Test group 9."""
    cases: list[dict] = []

    def _pairwise(ctrl: list[float], treated: list[tuple[int, list[float]]]) -> list[dict]:
        out = []
        ctrl_clean = [v for v in ctrl if v is not None and not (isinstance(v, float) and math.isnan(v))]
        if len(ctrl_clean) < 2:
            return []
        for dose, vals in treated:
            vc = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
            if len(vc) < 2:
                out.append({"dose_level": int(dose), "p_value_welch": None})
                continue
            _, p = sp_stats.ttest_ind(vc, ctrl_clean, equal_var=False)
            p_val = round(float(p), 6) if not math.isnan(p) else None
            out.append({"dose_level": int(dose), "p_value_welch": p_val})
        return out

    # 9.1 Rosetta-derived
    ctrl = ROSETTA_D1
    treated = [(10, ROSETTA_D2), (50, ROSETTA_D7[:5])]
    cases.append({
        "name": "two treated groups, Rosetta-derived",
        "category": "basic",
        "inputs": {"control": ctrl, "treated": [{"dose_level": d, "values": v} for d, v in treated]},
        "expected": _pairwise(ctrl, treated),
    })

    # 9.2 multiple treated groups
    ctrl2 = [2.0, 3.0, 4.0, 5.0]
    treated2 = [(10, [3.0, 4.0, 5.0, 6.0]), (50, [7.0, 8.0, 9.0, 10.0]), (100, [12.0, 13.0, 14.0, 15.0])]
    cases.append({
        "name": "three treated groups",
        "category": "basic",
        "inputs": {"control": ctrl2, "treated": [{"dose_level": d, "values": v} for d, v in treated2]},
        "expected": _pairwise(ctrl2, treated2),
    })

    # 9.3 raw p (not corrected) — same vals twice → same p
    ctrl3 = [1.0, 2.0, 3.0, 4.0]
    same = [3.0, 4.0, 5.0, 6.0]
    treated3 = [(10, same), (50, same)]
    cases.append({
        "name": "raw p-values (not multiplied by # tests)",
        "category": "basic",
        "inputs": {"control": ctrl3, "treated": [{"dose_level": d, "values": v} for d, v in treated3]},
        "expected": _pairwise(ctrl3, treated3),
    })

    # 9.4-9.5 NaN
    ctrl_nan = [None, 1.0, 2.0, 3.0]
    treated_nan = [(10, [10.0, 11.0, 12.0])]
    cases.append({
        "name": "NaN in control stripped",
        "category": "nan_handling",
        "inputs": {"control": ctrl_nan, "treated": [{"dose_level": d, "values": v} for d, v in treated_nan]},
        "expected": _pairwise(ctrl_nan, treated_nan),
    })

    treated_nan2 = [(10, [None, 10.0, 11.0, 12.0])]
    cases.append({
        "name": "NaN in treated stripped",
        "category": "nan_handling",
        "inputs": {"control": [1.0, 2.0, 3.0], "treated": [{"dose_level": d, "values": v} for d, v in treated_nan2]},
        "expected": _pairwise([1.0, 2.0, 3.0], treated_nan2),
    })

    # 9.6 control < 2 → []
    cases.append({
        "name": "control < 2 elements → empty list",
        "category": "edge_empty",
        "inputs": {"control": [1.0], "treated": [{"dose_level": 10, "values": [5.0, 6.0, 7.0]}]},
        "expected": [],
    })
    cases.append({
        "name": "no treated groups → empty list",
        "category": "edge_empty",
        "inputs": {"control": [1.0, 2.0, 3.0], "treated": []},
        "expected": [],
    })

    # 9.7 treated n=1 → None p
    cases.append({
        "name": "treated n=1 → p_value is None",
        "category": "edge_none",
        "inputs": {"control": [1.0, 2.0, 3.0], "treated": [{"dose_level": 10, "values": [5.0]}]},
        "expected": [{"dose_level": 10, "p_value_welch": None}],
    })

    _write("welch-pairwise", {"metadata": _meta("validate_rest_statistics.py — welch_pairwise"), "cases": cases})


def gen_bonferroni() -> None:
    """validate_rest_statistics.py — Test groups 10-11."""
    cases: list[dict] = []

    # Group 10: published
    cases.append({
        "name": "[Abdi] Table 1: 3 tests",
        "category": "published_reference",
        "inputs": {"p_values": ABDI_RAW_P, "n_tests": 3},
        "expected": [0.000120, 0.048300, 1.0],
    })
    cases.append({
        "name": "[Garcia/McDonald] 25-variable dietary",
        "category": "published_reference",
        "inputs": {"p_values": GARCIA_RAW_P, "n_tests": 25},
        "expected": [0.025, 0.200, 0.975, 1.0, 1.0, 1.0, 1.0],
    })

    # Group 11: algebraic
    cases.append({
        "name": "basic [0.01,0.04,0.05] × 3",
        "category": "algebraic",
        "inputs": {"p_values": [0.01, 0.04, 0.05]},
        "expected": [0.03, 0.12, 0.15],
    })
    cases.append({
        "name": "cap at 1.0",
        "category": "algebraic",
        "inputs": {"p_values": [0.5, 0.8]},
        "expected": [1.0, 1.0],
    })
    cases.append({
        "name": "explicit n_tests=5",
        "category": "algebraic",
        "inputs": {"p_values": [0.01, 0.02, 0.03], "n_tests": 5},
        "expected": [0.05, 0.10, 0.15],
    })
    cases.append({
        "name": "None passthrough (n_tests auto = 2)",
        "category": "edge_none_passthrough",
        "inputs": {"p_values": [0.01, None, 0.05]},
        "expected": [0.02, None, 0.10],
    })
    cases.append({
        "name": "all None",
        "category": "edge_none_passthrough",
        "inputs": {"p_values": [None, None, None]},
        "expected": [None, None, None],
    })
    cases.append({
        "name": "single p-value",
        "category": "edge",
        "inputs": {"p_values": [0.04]},
        "expected": [0.04],
    })
    cases.append({
        "name": "empty list",
        "category": "edge",
        "inputs": {"p_values": []},
        "expected": [],
    })
    cases.append({
        "name": "n_tests=0 returns originals",
        "category": "edge",
        "inputs": {"p_values": [0.01, 0.05], "n_tests": 0},
        "expected": [0.01, 0.05],
    })
    cases.append({
        "name": "order preserved",
        "category": "algebraic",
        "inputs": {"p_values": [0.05, 0.01, 0.03, 0.02]},
        "expected": [0.20, 0.04, 0.12, 0.08],
    })
    cases.append({
        "name": "very small p",
        "category": "algebraic",
        "inputs": {"p_values": [1e-10, 1e-8]},
        "expected": [2e-10, 2e-8],
    })

    _write("bonferroni", {"metadata": _meta("validate_rest_statistics.py — bonferroni_correct"), "cases": cases})


def gen_dunnett() -> None:
    """validate_dunnett.py — Dunnett 1955 dataset."""
    control = [7.40, 8.50, 7.20, 8.24, 9.84, 8.32]
    drug_a = [9.76, 8.80, 7.68, 9.36]
    drug_b = [12.80, 9.68, 12.16, 9.20, 10.55]

    ref = sp_stats.dunnett(np.array(drug_a), np.array(drug_b), control=np.array(control))
    expected = [
        {"dose_level": 1, "p_value": float(ref.pvalue[0]), "statistic": float(ref.statistic[0])},
        {"dose_level": 2, "p_value": float(ref.pvalue[1]), "statistic": float(ref.statistic[1])},
    ]

    cases = [{
        "name": "[Dunnett 1955] blood counts (millions cells/mm³)",
        "category": "published_dataset",
        "inputs": {
            "control": control,
            "treated": [{"dose_level": 1, "values": drug_a}, {"dose_level": 2, "values": drug_b}],
        },
        "expected": expected,
    }]

    # Add a few synthetic cross-checks
    syn_cases = [
        ("3 treated, balanced", [1.0, 2.0, 3.0, 2.5, 1.5], [(10, [3.0, 4.0, 5.0, 4.5]), (50, [5.0, 6.0, 7.0, 6.5]), (100, [7.0, 8.0, 9.0, 8.5])]),
        ("2 treated, unbalanced", [10.0, 11.0, 9.5, 10.5, 10.2, 11.3, 9.8], [(25, [12.0, 11.5, 13.0]), (75, [14.0, 15.0, 14.5, 13.5])]),
    ]
    for name, ctrl, treated in syn_cases:
        treated_arrays = [np.array(t[1]) for t in treated]
        ref = sp_stats.dunnett(*treated_arrays, control=np.array(ctrl))
        cases.append({
            "name": name,
            "category": "scipy_synthetic",
            "inputs": {
                "control": ctrl,
                "treated": [{"dose_level": d, "values": v} for d, v in treated],
            },
            "expected": [
                {"dose_level": int(treated[i][0]), "p_value": float(ref.pvalue[i]), "statistic": float(ref.statistic[i])}
                for i in range(len(treated))
            ],
        })

    _write("dunnett", {"metadata": _meta("validate_dunnett.py", {"p_value": TOL_PVALUE, "statistic": TOL_STATISTIC}), "cases": cases})


def gen_hedges_g() -> None:
    """validate_hedges_g.py — NIST Iris."""
    cases: list[dict] = []

    iris = [
        ("setosa (X=1)", IRIS_SETOSA_Y1, IRIS_SETOSA_Y2, 4.311260),
        ("versicolor (X=2)", IRIS_VERSI_Y1, IRIS_VERSI_Y2, 7.412040),
        ("virginica (X=3)", IRIS_VIRG_Y1, IRIS_VIRG_Y2, 7.168413),
    ]
    for name, y1, y2, nist_d in iris:
        a1, a2 = np.array(y1), np.array(y2)
        n1, n2 = len(a1), len(a2)
        df = n1 + n2 - 2
        pooled_std = math.sqrt(((n1 - 1) * float(np.var(a1, ddof=1)) + (n2 - 1) * float(np.var(a2, ddof=1))) / df)
        d_uncorrected = float((np.mean(a1) - np.mean(a2)) / pooled_std)
        j_approx = 1 - 3 / (4 * df - 1)
        j_exact = math.gamma(df / 2) / (math.sqrt(df / 2) * math.gamma((df - 1) / 2))
        g_with_correction = d_uncorrected * j_approx
        cases.append({
            "name": name,
            "category": "NIST_iris",
            "inputs": {"y1": y1, "y2": y2},
            "expected": {
                "cohens_d_uncorrected": d_uncorrected,
                "j_approx": j_approx,
                "j_exact": j_exact,
                "hedges_g_with_correction": g_with_correction,
                "df": df,
            },
            "published": {"nist_d_abs": nist_d, "tol": 1e-4},
        })

    # Edge cases
    cases.append({
        "name": "fewer than 2 in group → None",
        "category": "edge_none",
        "inputs": {"y1": [5.0], "y2": [1.0, 2.0, 3.0]},
        "expected": None,
    })
    cases.append({
        "name": "zero pooled std → None",
        "category": "edge_none",
        "inputs": {"y1": [5.0, 5.0, 5.0], "y2": [5.0, 5.0, 5.0]},
        "expected": None,
    })

    _write("hedges-g", {"metadata": _meta("validate_hedges_g.py"), "cases": cases})


def gen_jonckheere() -> None:
    """validate_trend_test.py — Jonckheere-Terpstra."""
    cases: list[dict] = []

    # Seeded synthetic data — we materialize the actual values
    rng = np.random.default_rng()  # not used; we use np.random for compat with original
    np.random.seed(42)
    g1 = np.random.normal(5, 1, size=10).tolist()
    g2 = np.random.normal(6, 1, size=12).tolist()
    g3 = np.random.normal(7, 1, size=8).tolist()
    groups = [g1, g2, g3]

    # Manual JT via pairwise Mann-Whitney U (the reference computation)
    J_ref = 0.0
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            u, _ = sp_stats.mannwhitneyu(groups[j], groups[i], alternative="two-sided")
            J_ref += float(u)
    N = sum(len(g) for g in groups)
    ns = [len(g) for g in groups]
    E_J = (N * N - sum(n * n for n in ns)) / 4.0
    Var_J = (N * N * (2 * N + 3) - sum(n * n * (2 * n + 3) for n in ns)) / 72.0
    Z_ref = float((J_ref - E_J) / math.sqrt(Var_J))
    p_ref = float(2 * (1 - sp_stats.norm.cdf(abs(Z_ref))))

    cases.append({
        "name": "[seed=42] 3 normal groups (5,6,7), n=(10,12,8)",
        "category": "seeded_synthetic",
        "data_origin": "np.random.seed(42); normal(5,1,10), normal(6,1,12), normal(7,1,8)",
        "inputs": {"groups": groups},
        "expected": {"statistic": Z_ref, "p_value": p_ref},
    })

    # Synthetic deterministic cases
    deterministic = [
        ("perfect monotone", [[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
        ("identical groups", [[1, 2, 3], [1, 2, 3], [1, 2, 3]]),
        ("reversed", [[7, 8, 9], [4, 5, 6], [1, 2, 3]]),
        ("two groups", [[1, 2, 3, 4], [5, 6, 7, 8]]),
    ]
    for name, groups_d in deterministic:
        # Compute manually since the JT is straightforward
        J = 0.0
        for i in range(len(groups_d)):
            for j in range(i + 1, len(groups_d)):
                arr_i = np.array(groups_d[i], dtype=float)
                arr_j = np.array(groups_d[j], dtype=float)
                diff = arr_j[:, None] - arr_i[None, :]
                J += float(np.sum(diff > 0)) + 0.5 * float(np.sum(diff == 0))
        N = sum(len(g) for g in groups_d)
        ns = [len(g) for g in groups_d]
        E_J = (N * N - sum(n * n for n in ns)) / 4.0
        Var_J = (N * N * (2 * N + 3) - sum(n * n * (2 * n + 3) for n in ns)) / 72.0
        if Var_J <= 0:
            cases.append({"name": name, "category": "deterministic", "inputs": {"groups": groups_d}, "expected": {"statistic": None, "p_value": None}})
            continue
        Z = (J - E_J) / math.sqrt(Var_J)
        p = 2 * (1 - sp_stats.norm.cdf(abs(Z)))
        cases.append({
            "name": name,
            "category": "deterministic",
            "inputs": {"groups": groups_d},
            "expected": {"statistic": float(Z), "p_value": float(p)},
        })

    # Edge cases
    cases.append({
        "name": "k < 2 → None",
        "category": "edge_none",
        "inputs": {"groups": [[1, 2, 3]]},
        "expected": {"statistic": None, "p_value": None},
    })
    cases.append({
        "name": "N < 4 → None",
        "category": "edge_none",
        "inputs": {"groups": [[1], [2]]},
        "expected": {"statistic": None, "p_value": None},
    })

    _write("jonckheere", {"metadata": _meta("validate_trend_test.py"), "cases": cases})


def gen_cochran_armitage() -> None:
    """validate_trend_test_incidence.py — basic CA test (binomial variance, default scores)."""
    cases: list[dict] = []

    # Group 5: Young Table 1 with default scores (0,1,…,k-1)
    young = [
        ("[Young] Table 1(a) 0/50, 0/50, 10/50", [0, 0, 10], [50, 50, 50]),
        ("[Young] Table 1(b) 0/50, 5/50, 10/50", [0, 5, 10], [50, 50, 50]),
        ("[Young] Table 1(c) 0/50, 10/50, 10/50", [0, 10, 10], [50, 50, 50]),
    ]
    for name, c, t in young:
        z, p = _ca_basic(c, t)
        cases.append({
            "name": name,
            "category": "young_table_1_default_scores",
            "inputs": {"counts": c, "totals": t},
            "expected": {"statistic": z, "p_value": p},
        })

    # Group 8: Tang rats
    z, p = _ca_basic(TANG_COUNTS, TANG_TOTALS)
    cases.append({
        "name": "[Tang 2006] 90-day neurotoxicity rats",
        "category": "published",
        "inputs": {"counts": TANG_COUNTS, "totals": TANG_TOTALS},
        "expected": {"statistic": z, "p_value": p},
    })

    # Group 9: DDT4 + Lilly Study 2
    z, p = _ca_basic(DDT4_COUNTS, DDT4_TOTALS)
    cases.append({
        "name": "[Young] Table 5c DDT 4-group",
        "category": "published",
        "inputs": {"counts": DDT4_COUNTS, "totals": DDT4_TOTALS},
        "expected": {"statistic": z, "p_value": p},
    })
    z, p = _ca_basic([10, 6, 7, 18], [120, 80, 80, 80])
    cases.append({
        "name": "[Young] Table 3b Lilly Study 2 males",
        "category": "published",
        "inputs": {"counts": [10, 6, 7, 18], "totals": [120, 80, 80, 80]},
        "expected": {"statistic": z, "p_value": p},
        "published": {"approx_z": 2.78, "tol": 0.01},
    })

    # Property tests / cross-checks (groups 6-7, 13)
    extra = [
        ("balanced 3-group", [5, 10, 20], [50, 50, 50]),
        ("balanced reversed", [20, 10, 5], [50, 50, 50]),
        ("flat 3-group", [10, 10, 10], [50, 50, 50]),
        ("4-group balanced", [10, 15, 20, 30], [100, 100, 100, 100]),
        ("4-group unbalanced", [3, 8, 15, 25], [40, 60, 80, 100]),
        ("strong trend large N", [10, 50, 100, 200], [1000, 1000, 1000, 1000]),
    ]
    for name, c, t in extra:
        z, p = _ca_basic(c, t)
        cases.append({
            "name": name,
            "category": "synthetic",
            "inputs": {"counts": c, "totals": t},
            "expected": {"statistic": z, "p_value": p},
        })

    # Group 14: flat large N
    flat_counts = [int(0.15 * 500)] * 4
    flat_totals = [500] * 4
    z, p = _ca_basic(flat_counts, flat_totals)
    cases.append({
        "name": "no trend, large N (flat 0.15 across 4 groups of 500)",
        "category": "synthetic",
        "inputs": {"counts": flat_counts, "totals": flat_totals},
        "expected": {"statistic": z, "p_value": p},
    })

    # Group 16: seeded random for k=2..6
    for k in [2, 3, 4, 5, 6]:
        np.random.seed(k)
        cc_k = np.random.randint(0, 20, k).tolist()
        cc_k[0] = 1
        cc_k[-1] = min(cc_k[-1], 49)
        tt_k = [50] * k
        z, p = _ca_basic(cc_k, tt_k)
        cases.append({
            "name": f"seeded k={k}",
            "category": "seeded_synthetic",
            "data_origin": f"np.random.seed({k}); randint(0,20,{k}); cc[0]=1, cc[-1]=min(cc[-1],49)",
            "inputs": {"counts": cc_k, "totals": tt_k},
            "expected": {"statistic": z, "p_value": p},
        })

    # Degenerate
    for name, c, t in [
        ("p̄ = 0", [0, 0, 0], [50, 50, 50]),
        ("p̄ = 1", [50, 50, 50], [50, 50, 50]),
        ("k = 1", [5], [50]),
        ("k = 0", [], []),
        ("N = 0", [0, 0], [0, 0]),
        ("single non-zero group (Sxx=0)", [0, 10], [0, 50]),
    ]:
        cases.append({
            "name": name,
            "category": "edge_none",
            "inputs": {"counts": c, "totals": t},
            "expected": {"statistic": None, "p_value": None},
        })

    _write("cochran-armitage", {"metadata": _meta("validate_trend_test_incidence.py"), "cases": cases})


def _ca_basic(counts: list[int], totals: list[int]) -> tuple[float | None, float | None]:
    """Reference implementation of basic CA trend test (binomial variance, scores=range(k))."""
    k = len(counts)
    if k < 2 or sum(totals) == 0:
        return None, None
    cc = np.array(counts, dtype=float)
    tt = np.array(totals, dtype=float)
    n = float(tt.sum())
    p_bar = cc.sum() / n
    if p_bar == 0 or p_bar == 1:
        return None, None
    d = np.arange(k, dtype=float)
    num = float(d @ cc - p_bar * (d @ tt))
    d_bar = float((d @ tt) / n)
    Sxx = float(sum(tt[i] * (d[i] - d_bar) ** 2 for i in range(k)))
    denom_sq = p_bar * (1 - p_bar) * Sxx
    if denom_sq <= 0:
        return None, None
    z = num / math.sqrt(denom_sq)
    p = 2 * (1 - float(sp_stats.norm.cdf(abs(z))))
    return float(z), float(p)


def gen_cochran_armitage_modified() -> None:
    """validate_trend_test_incidence_modified.py — extended CA + threshold_test."""
    cases_main: list[dict] = []
    cases_threshold: list[dict] = []

    # Group 3: Young Table 1 score sensitivity
    young_cases = [
        ("[Young] (a) [0,0,1]", YOUNG_A, [0, 0, 1]),
        ("[Young] (a) [0,1,2]", YOUNG_A, [0, 1, 2]),
        ("[Young] (a) [0,1,1]", YOUNG_A, [0, 1, 1]),
        ("[Young] (b) [0,0,1]", YOUNG_B, [0, 0, 1]),
        ("[Young] (b) [0,1,2]", YOUNG_B, [0, 1, 2]),
        ("[Young] (b) [0,1,1]", YOUNG_B, [0, 1, 1]),
        ("[Young] (c) [0,0,1]", YOUNG_C, [0, 0, 1]),
        ("[Young] (c) [0,1,2]", YOUNG_C, [0, 1, 2]),
        ("[Young] (c) [0,1,1]", YOUNG_C, [0, 1, 1]),
    ]
    for name, (c, t), sc in young_cases:
        z, p = _ca_full(c, t, scores=sc, alternative="increasing", variance="binomial")
        cases_main.append({
            "name": name,
            "category": "young_score_sensitivity",
            "inputs": {"counts": c, "totals": t, "scores": sc, "alternative": "increasing", "variance": "binomial"},
            "expected": {"z_statistic": z, "p_value": p},
        })

    # Group 4: DDT 7-group with 3 score sets
    ddt_cases = [
        ([0, 1, 1, 1, 1, 1, 1], 2.22),
        ([0, 1, 2, 3, 4, 5, 6], 5.25),
        ([0, 0, 0, 1, 2, 2, 2], 6.06),
    ]
    for sc, pub_z in ddt_cases:
        z, p = _ca_full(DDT_COUNTS, DDT_TOTALS, scores=sc, alternative="increasing", variance="binomial")
        cases_main.append({
            "name": f"[Young] Table 5b DDT scores={sc}",
            "category": "published_dose_response",
            "inputs": {"counts": DDT_COUNTS, "totals": DDT_TOTALS, "scores": sc, "alternative": "increasing", "variance": "binomial"},
            "expected": {"z_statistic": z, "p_value": p},
            "published": {"approx_z": pub_z, "tol": 0.20, "note": "1985 publication rounding"},
        })

    # Group 2: binomial vs hypergeometric on same data
    for label, c, t in [
        ("balanced 4-group", [10, 15, 20, 30], [100, 100, 100, 100]),
    ]:
        for var in ["binomial", "hypergeometric"]:
            z, p = _ca_full(c, t, scores=None, alternative="two-sided", variance=var)
            cases_main.append({
                "name": f"{label} variance={var}",
                "category": "variance_comparison",
                "inputs": {"counts": c, "totals": t, "alternative": "two-sided", "variance": var},
                "expected": {"z_statistic": z, "p_value": p},
            })

    # Group 7: affine invariance — same Z for different scores
    for sc in [[0, 1, 2], [1, 2, 3], [10, 30, 50]]:
        z, p = _ca_full([5, 10, 20], [50, 50, 50], scores=sc, alternative="two-sided", variance="binomial")
        cases_main.append({
            "name": f"affine invariance scores={sc}",
            "category": "affine_invariance",
            "inputs": {"counts": [5, 10, 20], "totals": [50, 50, 50], "scores": sc, "alternative": "two-sided", "variance": "binomial"},
            "expected": {"z_statistic": z, "p_value": p},
        })

    # Group 8: modified test
    for name, c, t in [
        ("uniform p̂", [10, 10, 10], [50, 50, 50]),
        ("heterogeneous p̂", [2, 10, 40], [100, 100, 100]),
        ("clear trend [5,15,30]", [5, 15, 30], [50, 50, 50]),
    ]:
        z, p = _ca_full(c, t, alternative="increasing", variance="binomial")
        z_mod, p_mod = _ca_modified(c, t, alternative="increasing")
        cases_main.append({
            "name": f"modified: {name}",
            "category": "modified_buonaccorsi",
            "inputs": {"counts": c, "totals": t, "alternative": "increasing", "modified": True},
            "expected": {"z_statistic": z, "p_value": p, "z_modified": z_mod, "p_value_modified": p_mod},
        })

    # Group 9: alternatives
    for alt in ["two-sided", "increasing", "decreasing"]:
        z, p = _ca_full([5, 10, 20], [50, 50, 50], alternative=alt, variance="binomial")
        cases_main.append({
            "name": f"alternative={alt}",
            "category": "alternative_directions",
            "inputs": {"counts": [5, 10, 20], "totals": [50, 50, 50], "alternative": alt, "variance": "binomial"},
            "expected": {"z_statistic": z, "p_value": p},
        })

    # Group 10: errors / degenerate
    cases_main.append({
        "name": "p̄ = 0 → degenerate (Z=0, p=1)",
        "category": "degenerate",
        "inputs": {"counts": [0, 0, 0], "totals": [50, 50, 50]},
        "expected": {"z_statistic": 0.0, "p_value": 1.0},
    })
    cases_main.append({
        "name": "p̄ = 1 → degenerate (Z=0, p=1)",
        "category": "degenerate",
        "inputs": {"counts": [50, 50, 50], "totals": [50, 50, 50]},
        "expected": {"z_statistic": 0.0, "p_value": 1.0},
    })
    cases_main.append({
        "name": "identical scores → degenerate",
        "category": "degenerate",
        "inputs": {"counts": [5, 10, 20], "totals": [50, 50, 50], "scores": [1, 1, 1]},
        "expected": {"z_statistic": 0.0, "p_value": 1.0},
    })
    for name, params in [
        ("k < 2", {"counts": [5], "totals": [50]}),
        ("count > total", {"counts": [60, 10], "totals": [50, 50]}),
        ("negative count", {"counts": [-1, 10], "totals": [50, 50]}),
        ("mismatched lengths", {"counts": [5, 10], "totals": [50, 50, 50]}),
        ("all totals 0", {"counts": [0, 0], "totals": [0, 0]}),
        ("invalid variance", {"counts": [5, 10], "totals": [50, 50], "variance": "invalid"}),
        ("invalid alternative", {"counts": [5, 10], "totals": [50, 50], "alternative": "wrong"}),
    ]:
        cases_main.append({
            "name": f"raises ValueError: {name}",
            "category": "validation_error",
            "inputs": params,
            "expected_error": "ValueError",
        })

    # Group 11: threshold_test
    ddt_4 = {"counts": [4, 1, 1, 7], "totals": [56, 49, 48, 41]}
    cases_threshold.append({
        "name": "[Young] Table 5c DDT 4-group, alpha=0.05, no Šidák",
        "category": "published",
        "inputs": {"counts": ddt_4["counts"], "totals": ddt_4["totals"], "alpha": 0.05, "adjust_alpha": False},
        "expected": _threshold_run(ddt_4["counts"], ddt_4["totals"], alpha=0.05, adjust_alpha=False),
        "published_zs": [-1.22, -0.79, 2.99],
        "published_effect_group": 3,
        "published_noel_groups": [0, 1, 2],
    })
    cases_threshold.append({
        "name": "[Young] DDT 4-group, alpha=0.05, with Šidák correction",
        "category": "published",
        "inputs": {"counts": ddt_4["counts"], "totals": ddt_4["totals"], "alpha": 0.05, "adjust_alpha": True},
        "expected": _threshold_run(ddt_4["counts"], ddt_4["totals"], alpha=0.05, adjust_alpha=True),
        "expected_alpha_adj": 1.0 - (1.0 - 0.05) ** (1.0 / 3.0),
    })
    cases_threshold.append({
        "name": "no trend → effect_group = None",
        "category": "no_effect",
        "inputs": {"counts": [5, 5, 5, 5], "totals": [50, 50, 50, 50], "alpha": 0.05, "adjust_alpha": False},
        "expected": _threshold_run([5, 5, 5, 5], [50, 50, 50, 50], alpha=0.05, adjust_alpha=False),
    })

    _write("cochran-armitage-modified", {
        "metadata": _meta("validate_trend_test_incidence_modified.py"),
        "cases": cases_main,
        "threshold_cases": cases_threshold,
    })


def _ca_full(counts: list[int], totals: list[int], scores: list[float] | None = None,
             alternative: str = "two-sided", variance: str = "binomial") -> tuple[float | None, float | None]:
    """Reference implementation of trend_test_incidence (modified, full signature)."""
    cc = np.array(counts, dtype=float)
    tt = np.array(totals, dtype=float)
    if cc.size == 0 or cc.size != tt.size:
        return None, None
    k = len(cc)
    if k < 2:
        return None, None
    n = float(tt.sum())
    if n == 0:
        return None, None
    p_bar = cc.sum() / n
    if p_bar == 0 or p_bar == 1:
        return 0.0, 1.0
    d = np.array(scores, dtype=float) if scores is not None else np.arange(k, dtype=float)
    num = float(d @ cc - p_bar * (d @ tt))
    d_bar = float((d @ tt) / n)
    Sxx = float(sum(tt[i] * (d[i] - d_bar) ** 2 for i in range(k)))
    if Sxx <= 0:
        return 0.0, 1.0
    if variance == "binomial":
        denom_sq = p_bar * (1 - p_bar) * Sxx
    else:
        denom_sq = p_bar * (1 - p_bar) * Sxx * n / (n - 1)
    z = num / math.sqrt(denom_sq)
    if alternative == "two-sided":
        p = 2 * float(sp_stats.norm.sf(abs(z)))
    elif alternative == "increasing":
        p = float(sp_stats.norm.sf(z))
    else:
        p = float(sp_stats.norm.cdf(z))
    return float(z), float(p)


def _ca_modified(counts: list[int], totals: list[int], alternative: str = "two-sided") -> tuple[float, float]:
    cc = np.array(counts, dtype=float)
    tt = np.array(totals, dtype=float)
    k = len(cc)
    n = float(tt.sum())
    p_bar = cc.sum() / n
    d = np.arange(k, dtype=float)
    num = float(d @ cc - p_bar * (d @ tt))
    d_bar = float((d @ tt) / n)
    devs = d - d_bar
    p_hat = np.divide(cc, tt, out=np.zeros_like(cc), where=tt > 0)
    sigma2_m = float(np.sum(tt * devs ** 2 * p_hat * (1.0 - p_hat)))
    if sigma2_m <= 0:
        return 0.0, 1.0
    z = num / math.sqrt(sigma2_m)
    if alternative == "two-sided":
        p = 2 * float(sp_stats.norm.sf(abs(z)))
    elif alternative == "increasing":
        p = float(sp_stats.norm.sf(z))
    else:
        p = float(sp_stats.norm.cdf(z))
    return float(z), float(p)


def _threshold_run(counts: list[int], totals: list[int], alpha: float = 0.05, adjust_alpha: bool = True) -> list[dict]:
    """Reference implementation of threshold_test (Williams-style sequential)."""
    counts = list(counts)
    totals = list(totals)
    k = len(counts)
    n_comparisons = k - 1
    alpha_adj = 1.0 - (1.0 - alpha) ** (1.0 / n_comparisons) if adjust_alpha else alpha

    results: list[dict] = []
    pool_count, pool_total = counts[0], totals[0]
    for i in range(1, k):
        z, p = _ca_full([pool_count, counts[i]], [pool_total, totals[i]], scores=[0, 1], alternative="increasing")
        sig = p is not None and p <= alpha_adj
        entry = {
            "test": f"groups {list(range(i))} vs group {i}",
            "control_count": pool_count,
            "control_total": pool_total,
            "control_pct": pool_count / pool_total * 100 if pool_total else 0,
            "treated_count": counts[i],
            "treated_total": totals[i],
            "treated_pct": counts[i] / totals[i] * 100 if totals[i] else 0,
            "z": z,
            "p": p,
            "alpha_adj": alpha_adj,
            "significant": sig,
        }
        if sig:
            entry["effect_group"] = i
            entry["noel_groups"] = list(range(i))
            results.append(entry)
            return results
        pool_count += counts[i]
        pool_total += totals[i]
        results.append(entry)
    if results:
        results[-1]["noel_groups"] = list(range(k))
        results[-1]["effect_group"] = None
    return results


# ── Fisher exact (independent oracle, no scipy in oracle) ──

def _hypergeom_pmf(a: int, N: int, R0: int, C0: int) -> float:
    R1 = N - R0
    b = R0 - a
    c = C0 - a
    d = R1 - c
    if any(x < 0 for x in (a, b, c, d)):
        return 0.0
    return comb(R0, a) * comb(R1, c) / comb(N, C0)


def _fisher_p_two_sided(table: list[list[int]]) -> float:
    a, b = table[0]
    c, d = table[1]
    R0, R1 = a + b, c + d
    C0 = a + c
    N = R0 + R1
    a_min = max(0, C0 - R1)
    a_max = min(R0, C0)
    p_obs = _hypergeom_pmf(a, N, R0, C0)
    return sum(_hypergeom_pmf(x, N, R0, C0) for x in range(a_min, a_max + 1) if _hypergeom_pmf(x, N, R0, C0) <= p_obs + 1e-15)


def _fisher_p_greater(table: list[list[int]]) -> float:
    a, b = table[0]
    c, d = table[1]
    R0, R1 = a + b, c + d
    C0 = a + c
    N = R0 + R1
    a_max = min(R0, C0)
    return sum(_hypergeom_pmf(x, N, R0, C0) for x in range(a, a_max + 1))


def _fisher_p_less(table: list[list[int]]) -> float:
    a, b = table[0]
    c, d = table[1]
    R0, R1 = a + b, c + d
    C0 = a + c
    N = R0 + R1
    a_min = max(0, C0 - R1)
    return sum(_hypergeom_pmf(x, N, R0, C0) for x in range(a_min, a + 1))


def _odds_ratio(table: list[list[int]]) -> float:
    a, b = table[0]
    c, d = table[1]
    if b * c == 0:
        return float("inf") if a * d > 0 else 0.0
    return (a * d) / (b * c)


def gen_fisher_exact() -> None:
    """validate_fisher_boschloo.py — Fisher portion (skip Boschloo: not in our port)."""
    cases: list[dict] = []

    # Reference tables with hand-computed exact p-values
    hand_tables = {
        "[Fisher35] Lady Tea [[3,1],[1,3]]": ([[3, 1], [1, 3]], {"p_two_exact": 34/70, "p_gt_exact": 17/70, "or_expected": 9.0}),
        "[Agresti02] [[6,2],[1,4]]": ([[6, 2], [1, 4]], {"p_two_exact": 176/1716, "p_gt_exact": 148/1716, "or_expected": 12.0}),
        "Toxicology [[8,2],[2,8]]": ([[8, 2], [2, 8]], {"or_expected": 16.0}),
        "Null [[5,5],[5,5]]": ([[5, 5], [5, 5]], {"or_expected": 1.0}),
        "[Saari04] [[74,31],[43,32]]": ([[74, 31], [43, 32]], {"or_expected": (74 * 32) / (31 * 43)}),
        "[Rosner06] [[3,17],[1,19]]": ([[3, 17], [1, 19]], {"or_expected": (3 * 19) / (17 * 1)}),
        "Perfect [[5,0],[0,5]]": ([[5, 0], [0, 5]], {"p_two_exact": 1/126}),
        "Inverse [[0,5],[5,0]]": ([[0, 5], [5, 0]], {}),
        "Extreme [[10,0],[0,10]]": ([[10, 0], [0, 10]], {"p_two_exact": 1/92378}),
        "Weak [[4,3],[2,5]]": ([[4, 3], [2, 5]], {}),
        "Large balanced [[30,10],[10,30]]": ([[30, 10], [10, 30]], {}),
        "Unequal margins [[1,9],[5,5]]": ([[1, 9], [5, 5]], {}),
    }
    for name, (tbl, extra) in hand_tables.items():
        oracle_p_two = _fisher_p_two_sided(tbl)
        oracle_p_gt = _fisher_p_greater(tbl)
        oracle_p_lt = _fisher_p_less(tbl)
        oracle_or = _odds_ratio(tbl)
        case: dict = {
            "name": name,
            "category": "published_or_constructed",
            "inputs": {"table": tbl},
            "expected": {
                "p_two_sided": oracle_p_two,
                "p_greater": oracle_p_gt,
                "p_less": oracle_p_lt,
                "odds_ratio": oracle_or if oracle_or != float("inf") else "Infinity",
            },
        }
        if extra:
            case["published"] = extra
        cases.append(case)

    # Group 7: 100 random tables (seeded)
    np.random.seed(42)
    n_cross = 100
    random_cases = []
    for i in range(n_cross):
        tbl = np.random.randint(0, 25, size=(2, 2)).tolist()
        # Skip degenerate (all-zero row or col)
        if sum(tbl[0]) == 0 or sum(tbl[1]) == 0:
            continue
        if tbl[0][0] + tbl[1][0] == 0 or tbl[0][1] + tbl[1][1] == 0:
            continue
        oracle_p_two = _fisher_p_two_sided(tbl)
        random_cases.append({
            "table": tbl,
            "expected_p_two_sided": oracle_p_two,
        })

    cases.append({
        "name": f"random {len(random_cases)} valid 2×2 tables (seeded)",
        "category": "randomized_property",
        "data_origin": f"np.random.seed(42); randint(0,25,(2,2)) × {n_cross}, dropping degenerate",
        "tables": random_cases,
    })

    _write("fisher-exact", {"metadata": _meta("validate_fisher_boschloo.py — Fisher portion"), "cases": cases})


def gen_williams() -> None:
    """validate_fixed_williams.py — uses published 1971/1972 paper data."""
    cases: list[dict] = []

    # Group 1: PAVA tests (paper data)
    cases.append({
        "name": "[1971] §2 p.106 PAVA increasing on all 7 groups",
        "category": "pava",
        "inputs": {
            "values": [10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7],
            "weights": [1, 1, 1, 1, 1, 1, 1],
            "direction": "increase",
        },
        "expected": [10.1, 10.1, 10.1, 10.6, 11.4, 11.8, 11.8],
        "tolerance": 0.05,
    })
    cases.append({
        "name": "[1972] §7 p.530 PAVA decreasing on dose groups",
        "category": "pava",
        "inputs": {
            "values": [52.7, 45.2, 47.1, 44.8, 46.6],
            "weights": [10, 9, 10, 10, 8],
            "direction": "decrease",
        },
        "expected": [52.7, 46.2, 46.2, 45.6, 45.6],
        "tolerance": 0.05,
    })
    # Edge PAVA cases
    cases.append({
        "name": "PAVA monotone input → unchanged",
        "category": "pava",
        "inputs": {"values": [1.0, 2.0, 3.0, 4.0], "weights": [1, 1, 1, 1], "direction": "increase"},
        "expected": [1.0, 2.0, 3.0, 4.0],
    })
    cases.append({
        "name": "PAVA constant input → unchanged",
        "category": "pava",
        "inputs": {"values": [5.0, 5.0, 5.0], "weights": [1, 1, 1], "direction": "increase"},
        "expected": [5.0, 5.0, 5.0],
    })
    cases.append({
        "name": "PAVA fully reversed → global mean",
        "category": "pava",
        "inputs": {"values": [4.0, 3.0, 2.0, 1.0], "weights": [1, 1, 1, 1], "direction": "increase"},
        "expected": [2.5, 2.5, 2.5, 2.5],
    })

    # Group 3: 1971 critical-value spot checks
    cv_1971 = [
        ("k=6, v=42, α=0.05", {"k": 6, "v": 42, "alpha": 0.05}, 1.81),
        ("k=6, v=42, α=0.01", {"k": 6, "v": 42, "alpha": 0.01}, 2.50),
        ("k=5, v=42, α=0.05", {"k": 5, "v": 42, "alpha": 0.05}, 1.80),
        ("k=4, v=42, α=0.05", {"k": 4, "v": 42, "alpha": 0.05}, 1.80),
        ("k=3, v=42, α=0.05", {"k": 3, "v": 42, "alpha": 0.05}, 1.79),
        ("k=1, v=∞, α=0.05", {"k": 1, "v": "inf", "alpha": 0.05}, 1.645),
        ("k=10, v=5, α=0.05", {"k": 10, "v": 5, "alpha": 0.05}, 2.25),
        ("k=1, v=5, α=0.05", {"k": 1, "v": 5, "alpha": 0.05}, 2.02),
        ("k=1, v=∞, α=0.01", {"k": 1, "v": "inf", "alpha": 0.01}, 2.326),
        ("k=2, v=5, α=0.01", {"k": 2, "v": 5, "alpha": 0.01}, 3.50),
    ]
    for name, params, expected in cv_1971:
        cases.append({
            "name": f"[1971] {name}",
            "category": "lookup_1971",
            "inputs": params,
            "expected": expected,
        })

    # Group 4: 1972 critical-value spot checks
    cv_1972 = [
        ("i=2, v=60, α=0.025, w=2", {"i": 2, "v": 60, "alpha": 0.025, "w": 18/9}, 2.05, 0.02),
        ("i=3, v=60, α=0.025, w=1.8", {"i": 3, "v": 60, "alpha": 0.025, "w": 18/10}, 2.06, 0.02),
        ("i=4, v=60, α=0.025, w=1.8", {"i": 4, "v": 60, "alpha": 0.025, "w": 18/10}, 2.07, 0.02),
        ("i=5, v=60, α=0.025, w=2.25", {"i": 5, "v": 60, "alpha": 0.025, "w": 18/8}, 2.06, 0.02),
        ("i=6, v=∞, α=0.050, w=1", {"i": 6, "v": "inf", "alpha": 0.050, "w": 1.0}, 1.760, 0.0001),
        ("i=2, v=5, α=0.010, w=1", {"i": 2, "v": 5, "alpha": 0.010, "w": 1.0}, 3.501, 0.0001),
        ("i=10, v=∞, α=0.005, w=1", {"i": 10, "v": "inf", "alpha": 0.005, "w": 1.0}, 2.623, 0.0001),
    ]
    for name, params, expected, tol in cv_1972:
        cases.append({
            "name": f"[1972] {name}",
            "category": "lookup_1972",
            "inputs": params,
            "expected": expected,
            "tolerance": tol,
        })

    # Group 5: full integration — paper examples
    # 1971 §2+§5 example
    cases.append({
        "name": "[1971] §5 p.117 full step-down (k=6, r=8, s²=1.16)",
        "category": "integration",
        "inputs": {
            "means": [10.4, 9.9, 10.0, 10.6, 11.4, 11.9, 11.7],
            "sds": [math.sqrt(1.16)] * 7,
            "ns": [8, 8, 8, 8, 8, 8, 8],
            "dose_labels": ["0", "1", "2", "3", "4", "5", "6"],
            "direction": "increase",
            "alpha": 0.05,
        },
        "expected": {
            "minimum_effective_dose": "4",
            "n_step_down_results": 4,
            "first_three_significant": True,
            "fourth_significant": False,
            "t_dose6": 2.60,
            "t_dose3": 0.37,
            "tol": 0.01,
        },
    })

    # 1972 §7 example
    cases.append({
        "name": "[1972] §7 p.530 full step-down (k=5, decrease, α=0.025)",
        "category": "integration",
        "inputs": {
            "means": [50.1, 52.7, 45.2, 47.1, 44.8, 46.6],
            "sds": [math.sqrt(22.28)] * 6,
            "ns": [18, 10, 9, 10, 10, 8],
            "dose_labels": ["0", "1", "2", "3", "4", "5"],
            "direction": "decrease",
            "alpha": 0.025,
        },
        "expected": {
            "constrained_means": [50.1, 52.7, 46.2, 46.2, 45.6, 45.6],
            "minimum_effective_dose": "3",
            "stats_by_dose_index": {"5": 2.24, "4": 2.42, "3": 2.09, "2": 2.02},
            "dose2_significant": False,
            "tol_means": 0.05,
            "tol_stats": 0.02,
        },
    })

    # Group 6 edge cases
    cases.append({
        "name": "k=1 reduces to Student's t",
        "category": "integration_edge",
        "inputs": {"means": [10.0, 12.0], "sds": [1.0, 1.0], "ns": [10, 10], "dose_labels": ["0", "1"], "direction": "increase", "alpha": 0.05},
        "expected": {"critical_value_step1": float(sp_stats.t.ppf(0.95, 18)), "tol_cv": 0.001},
    })
    cases.append({
        "name": "no effect: all means equal → MED = None",
        "category": "integration_edge",
        "inputs": {"means": [10.0, 10.0, 10.0, 10.0], "sds": [1.0, 1.0, 1.0, 1.0], "ns": [10, 10, 10, 10], "dose_labels": ["0", "1", "2", "3"], "direction": "increase", "alpha": 0.05},
        "expected": {"minimum_effective_dose": None},
    })
    cases.append({
        "name": "strong effect: all doses significant → MED = 1",
        "category": "integration_edge",
        "inputs": {"means": [0.0, 5.0, 10.0, 15.0], "sds": [1.0, 1.0, 1.0, 1.0], "ns": [20, 20, 20, 20], "dose_labels": ["0", "1", "2", "3"], "direction": "increase", "alpha": 0.05},
        "expected": {"minimum_effective_dose": "1", "all_groups_tested": True},
    })
    cases.append({
        "name": "auto-detect direction: highest < control → 'decrease'",
        "category": "integration_edge",
        "inputs": {"means": [50.0, 45.0, 40.0], "sds": [2.0, 2.0, 2.0], "ns": [10, 10, 10], "dose_labels": ["0", "1", "2"], "direction": "auto", "alpha": 0.05},
        "expected": {"direction": "decrease"},
    })
    cases.append({
        "name": "auto-detect direction: highest > control → 'increase'",
        "category": "integration_edge",
        "inputs": {"means": [10.0, 12.0, 15.0], "sds": [2.0, 2.0, 2.0], "ns": [10, 10, 10], "dose_labels": ["0", "1", "2"], "direction": "auto", "alpha": 0.05},
        "expected": {"direction": "increase"},
    })

    _write("williams", {"metadata": _meta("validate_fixed_williams.py"), "cases": cases})


def gen_ancova() -> None:
    """validate_ancova.py — Montgomery 15.10 + ESS apple trees."""
    cases: list[dict] = []

    # Run scipy-equivalent OLS on Montgomery to capture all expected outputs
    montgomery_result = _ancova_run(MONTGOMERY_STR, MONTGOMERY_DIA, MONTGOMERY_GRP, control_group=1)
    cases.append({
        "name": "[Montgomery 15.10] / [Purdue/SAS] — fiber strength by machine, diameter as covariate",
        "category": "published_dataset",
        "inputs": {
            "organ_values": MONTGOMERY_STR,
            "body_weights": MONTGOMERY_DIA,
            "groups": MONTGOMERY_GRP,
            "control_group": 1,
            "use_organ_free_bw": False,
            "alpha": 0.05,
        },
        "expected": montgomery_result,
        "published": {
            "sas_r_squared": 0.919209,
            "sas_mse": 2.5441718,
            "sas_ls_means": {"1": 40.3824131, "2": 41.4192229, "3": 38.7983640},
            "sas_slope_t": 8.365,
            "sas_homogeneity_f": 0.49,
            "sas_homogeneity_p": 0.6293,
            "df_error": 11,
            "df_interaction": 9,
            "tol_r2": 0.001,
            "tol_mse": 0.001,
            "tol_ls_mean": 0.001,
            "tol_slope_t": 0.01,
            "tol_homog_f": 0.01,
            "tol_homog_p": 0.001,
        },
    })

    # ESS apple trees (CRD)
    ess_result = _ancova_run(ESS_Y, ESS_X, ESS_GRP, control_group=0)
    cases.append({
        "name": "[ESS] p.128 — apple trees CRD",
        "category": "published_dataset",
        "inputs": {
            "organ_values": ESS_Y,
            "body_weights": ESS_X,
            "groups": ESS_GRP,
            "control_group": 0,
            "use_organ_free_bw": False,
            "alpha": 0.05,
        },
        "expected": ess_result,
        "published": {
            "x_overall_mean": 8.308,
            "tol_x_mean": 0.001,
        },
    })

    # Lazic organ-free covariate
    lazic_result = _ancova_run(MONTGOMERY_STR, MONTGOMERY_DIA, MONTGOMERY_GRP, control_group=1, use_organ_free_bw=True)
    cases.append({
        "name": "[Lazic 2020] organ-free body weight covariate",
        "category": "lazic_variant",
        "inputs": {
            "organ_values": MONTGOMERY_STR,
            "body_weights": MONTGOMERY_DIA,
            "groups": MONTGOMERY_GRP,
            "control_group": 1,
            "use_organ_free_bw": True,
            "alpha": 0.05,
        },
        "expected": lazic_result,
        "expected_covariate_mean": float(np.mean(np.array(MONTGOMERY_DIA, dtype=float) - np.array(MONTGOMERY_STR, dtype=float))),
    })

    # Edge cases
    cases.append({
        "name": "insufficient data (n=2, k=2)",
        "category": "edge_none",
        "inputs": {"organ_values": [1.0, 2.0], "body_weights": [3.0, 4.0], "groups": [0, 1], "control_group": 0},
        "expected": None,
    })
    cases.append({
        "name": "single group → None",
        "category": "edge_none",
        "inputs": {"organ_values": [1.0, 2.0, 3.0], "body_weights": [4.0, 5.0, 6.0], "groups": [0, 0, 0], "control_group": 0},
        "expected": None,
    })

    # Seeded equal means → no significant pairwise differences
    np.random.seed(99)
    y_eq = np.random.normal(100, 5, 30).tolist()
    x_eq = np.random.normal(50, 10, 30).tolist()
    g_eq = ([0] * 10) + ([1] * 10) + ([2] * 10)
    eq_result = _ancova_run(y_eq, x_eq, g_eq, control_group=0)
    cases.append({
        "name": "seeded equal means (np.random.seed=99, n=30) → no significant differences",
        "category": "seeded_property",
        "data_origin": "np.random.seed(99); normal(100,5,30); normal(50,10,30); groups [0]*10+[1]*10+[2]*10",
        "inputs": {
            "organ_values": y_eq,
            "body_weights": x_eq,
            "groups": g_eq,
            "control_group": 0,
            "alpha": 0.05,
        },
        "expected": eq_result,
        "expected_property": "all_pairwise_p_above_0.05",
    })

    _write("ancova", {"metadata": _meta("validate_ancova.py"), "cases": cases})


def _ancova_run(organ_values: list[float], body_weights: list[float], groups: list[int],
                control_group: int = 0, use_organ_free_bw: bool = False, alpha: float = 0.05) -> dict | None:
    """Reference implementation of run_ancova matching ancova.py exactly."""
    ov = np.array(organ_values, dtype=float)
    bw = np.array(body_weights, dtype=float)
    gp = np.array(groups, dtype=int)
    mask = ~(np.isnan(ov) | np.isnan(bw))
    ov, bw, gp = ov[mask], bw[mask], gp[mask]

    unique_groups = sorted(set(gp.tolist()))
    k = len(unique_groups)
    n = len(ov)
    if n < k + 2 or k < 2:
        return None

    cov = (bw - ov) if use_organ_free_bw else bw.copy()
    cov_mean = float(np.mean(cov))

    treated = [g for g in unique_groups if g != control_group]
    p_a = 1 + len(treated) + 1
    X_a = np.zeros((n, p_a))
    X_a[:, 0] = 1
    for j, g in enumerate(treated):
        X_a[:, 1 + j] = (gp == g).astype(float)
    X_a[:, -1] = cov

    beta_a, _, _, _ = np.linalg.lstsq(X_a, ov, rcond=None)
    resid = ov - X_a @ beta_a
    rss_a = float(np.sum(resid ** 2))
    df_a = n - p_a
    mse = rss_a / df_a if df_a > 0 else float("inf")
    try:
        XtX_inv = np.linalg.inv(X_a.T @ X_a)
    except np.linalg.LinAlgError:
        XtX_inv = np.linalg.pinv(X_a.T @ X_a)
    vcov_a = mse * XtX_inv

    p_int = p_a + len(treated)
    X_int = np.zeros((n, p_int))
    X_int[:, :p_a] = X_a
    for j, g in enumerate(treated):
        X_int[:, p_a + j] = X_a[:, 1 + j] * cov
    beta_int, _, _, _ = np.linalg.lstsq(X_int, ov, rcond=None)
    resid_int = ov - X_int @ beta_int
    rss_int = float(np.sum(resid_int ** 2))
    df_int = n - p_int

    if df_a - df_int > 0 and df_int > 0 and rss_int > 0:
        f_hom = ((rss_a - rss_int) / (df_a - df_int)) / (rss_int / df_int)
        p_hom = float(1 - sp_stats.f.cdf(f_hom, df_a - df_int, df_int))
    else:
        f_hom, p_hom = None, None
    homogeneous = p_hom is None or p_hom >= alpha

    slope_est = float(beta_a[-1])
    slope_se = float(math.sqrt(max(0.0, float(vcov_a[-1, -1]))))
    slope_t = slope_est / slope_se if slope_se > 0 else 0.0
    slope_p = float(2 * (1 - sp_stats.t.cdf(abs(slope_t), df_a)))

    tss = float(np.sum((ov - np.mean(ov)) ** 2))
    r_sq = 1 - rss_a / tss if tss > 0 else 0.0

    adjusted_means = []
    raw_means: dict[int, float] = {}
    for g in unique_groups:
        g_mask = gp == g
        g_vals = ov[g_mask]
        raw_mean = float(np.mean(g_vals))
        raw_means[g] = raw_mean
        g_n = int(np.sum(g_mask))
        x_pred = np.zeros(p_a)
        x_pred[0] = 1
        if g in treated:
            x_pred[1 + treated.index(g)] = 1
        x_pred[-1] = cov_mean
        adj_mean = float(x_pred @ beta_a)
        adj_se = float(math.sqrt(max(0.0, float(x_pred @ vcov_a @ x_pred))))
        adjusted_means.append({
            "group": int(g),
            "raw_mean": round(raw_mean, 4),
            "adjusted_mean": round(adj_mean, 4),
            "n": g_n,
            "se": round(adj_se, 4),
        })

    pairwise = []
    for g in treated:
        j = treated.index(g)
        c_vec = np.zeros(p_a)
        c_vec[1 + j] = 1
        diff = float(c_vec @ beta_a)
        se_diff = float(math.sqrt(max(0.0, float(c_vec @ vcov_a @ c_vec))))
        t_stat = diff / se_diff if se_diff > 0 else 0.0
        p_val = float(2 * (1 - sp_stats.t.cdf(abs(t_stat), df_a)))
        pairwise.append({
            "group": int(g),
            "difference": round(diff, 4),
            "se": round(se_diff, 4),
            "t_statistic": round(t_stat, 4),
            "p_value": round(p_val, 6),
            "significant": p_val < alpha,
        })

    sqrt_mse = math.sqrt(mse) if mse > 0 else 1.0
    ctrl_raw = raw_means.get(control_group, 0.0)
    effect_decomp = []
    for pw in pairwise:
        g = pw["group"]
        total = raw_means[g] - ctrl_raw
        direct = pw["difference"]
        indirect = total - direct
        prop_direct = direct / total if abs(total) > 1e-10 else 1.0
        direct_d = direct / sqrt_mse if sqrt_mse > 0 else 0.0
        j_corr = 1 - 3 / (4 * df_a - 1) if df_a > 1 else 1.0
        direct_g = abs(direct_d * j_corr)
        effect_decomp.append({
            "group": int(g),
            "total_effect": round(total, 4),
            "direct_effect": round(direct, 4),
            "indirect_effect": round(indirect, 4),
            "proportion_direct": round(prop_direct, 4),
            "direct_g": round(direct_g, 4),
            "direct_p": pw["p_value"],
        })

    return {
        "adjusted_means": adjusted_means,
        "pairwise": pairwise,
        "slope": {
            "estimate": round(slope_est, 6),
            "se": round(slope_se, 6),
            "t_statistic": round(slope_t, 4),
            "p_value": round(slope_p, 6),
        },
        "slope_homogeneity": {
            "f_statistic": round(f_hom, 4) if f_hom is not None else None,
            "p_value": round(p_hom, 6) if p_hom is not None else None,
            "homogeneous": homogeneous,
        },
        "effect_decomposition": effect_decomp,
        "model_r_squared": round(r_sq, 4),
        "mse": round(mse, 6),
        "use_organ_free_bw": use_organ_free_bw,
        "covariate_mean": round(cov_mean, 4),
    }


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════

def main() -> None:
    print(f"Generating fixtures into {FIXTURES_DIR}")
    print()
    gen_welch_t()
    gen_mann_whitney()
    gen_spearman()
    gen_severity_trend()
    gen_welch_pairwise()
    gen_bonferroni()
    gen_dunnett()
    gen_hedges_g()
    gen_jonckheere()
    gen_cochran_armitage()
    gen_cochran_armitage_modified()
    gen_fisher_exact()
    gen_williams()
    gen_ancova()
    print()
    print("Done.")


if __name__ == "__main__":
    main()
