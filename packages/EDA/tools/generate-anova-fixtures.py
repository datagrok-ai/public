"""
Generate anova-fixtures.ts from NIST StRD datasets and scipy reference values.

This script is OFFLINE TOOLING — it is not part of the package bundle.
Run it whenever the fixtures need to be regenerated:

    python generate-anova-fixtures.py

Requirements:
    scipy >= 1.10  (works with 1.16, 1.17)
    numpy
    Network access to https://www.itl.nist.gov (with browser User-Agent)

Output:
    ../src/tests/anova-fixtures.ts

Datasets (NIST StRD, https://www.itl.nist.gov/div898/strd/anova/):
    SiRstv  — Lower difficulty,  25 obs,   5 groups
    AtmWtAg — Average difficulty, 48 obs,  2 groups
    SmLs04  — Average difficulty, 1800 obs, 9 groups (5 const digits)
    SmLs07  — Higher difficulty,  18 obs,   9 groups (11 const digits)
    SmLs09  — Higher difficulty,  18009 obs, 9 groups (11 const digits)

Plus a small validation_simple dataset (the 15-point one mirroring the
existing TS test fixture).
"""

import json
import os
import ssl
import sys
import urllib.request
from pathlib import Path

import numpy as np
from scipy import stats


NIST_BASE = "https://www.itl.nist.gov/div898/strd/anova"
DATASETS = ["SiRstv", "AtmWtAg", "SmLs04", "SmLs07", "SmLs09"]

# Per-dataset tolerances {rtol, atol}.
#
# Fisher F/SS reference comes from scipy (which itself matches NIST certified
# values up to scipy's own precision); Welch reference always comes from scipy.
# Welford accumulation in Float64 matches scipy near ULP for normal data
# (SiRstv, AtmWtAg, SmLs04) but loses meaningful precision on the deliberately
# stiff NIST datasets (SmLs07: 11 const digits; SmLs09: 11 const digits × 2001
# obs/group). The looser tolerances on those two record observed accuracy and
# serve as a regression guard against further degradation.
TOLERANCE = {
    "SiRstv":  {"fisher": {"rtol": 1e-10, "atol": 1e-12}, "welch": {"rtol": 1e-9, "atol": 1e-12}},
    "AtmWtAg": {"fisher": {"rtol": 1e-7,  "atol": 1e-9 }, "welch": {"rtol": 1e-7, "atol": 1e-9}},
    "SmLs04":  {"fisher": {"rtol": 1e-7,  "atol": 1e-9 }, "welch": {"rtol": 1e-7, "atol": 1e-9}},
    "SmLs07":  {"fisher": {"rtol": 1e-2,  "atol": 1e-3 }, "welch": {"rtol": 1e-2, "atol": 1e-3}},
    "SmLs09":  {"fisher": {"rtol": 1e-3,  "atol": 1e-3 }, "welch": {"rtol": 1e-4, "atol": 1e-3}},
}
SIMPLE_TOLERANCE = {"rtol": 1e-6, "atol": 1e-9}


def _ssl_ctx():
    ctx = ssl.create_default_context()
    # NIST cert has a basic-constraint quirk; bypass under our corp env.
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def fetch(name: str) -> str:
    req = urllib.request.Request(
        f"{NIST_BASE}/{name}.dat",
        headers={"User-Agent": "Mozilla/5.0"},
    )
    return urllib.request.urlopen(req, timeout=30, context=_ssl_ctx()).read().decode("ascii")


def parse_nist_dat(text: str):
    """Parse a NIST StRD ANOVA .dat file.

    Returns dict with keys:
        factor:  list[int]    — 1-based group label
        value:   list[float]  — observation
        cert_f:  float        — certified F-statistic (or None)
        cert_ssb: float       — certified between-group SS (or None)
        cert_ssw: float       — certified within-group SS (or None)
        data_start_line: int  — 1-based line where data starts
    """
    lines = text.splitlines()

    # Find "Data" line range from header.
    data_start = None
    for line in lines[:60]:
        if "Data" in line and "lines" in line and "to" in line:
            # e.g.  "                Data               (lines 61 to 85)"
            try:
                inside = line[line.index("("):]
                nums = [int(t) for t in inside.replace("(", " ").replace(")", " ")
                        .replace(",", " ").split() if t.isdigit()]
                if len(nums) >= 2:
                    data_start = nums[0]
                    break
            except Exception:
                pass
    if data_start is None:
        raise RuntimeError("could not locate data line range in header")

    # Pull certified F, SSB, SSW out of the table.
    cert_f, cert_ssb, cert_ssw = None, None, None
    for line in lines:
        s = line.strip()
        if s.startswith("Between"):
            # "Between Instrument  4  5.114...E-02  1.27865...E-02  1.18046...E+00"
            parts = s.split()
            # last is F, second-to-last is MS, third-to-last is SS, before is df
            try:
                cert_f = float(parts[-1])
                cert_ssb = float(parts[-3])
            except Exception:
                pass
        elif s.startswith("Within"):
            parts = s.split()
            try:
                cert_ssw = float(parts[-2])
            except Exception:
                pass

    # Parse data rows.
    factor: list[int] = []
    value: list[float] = []
    for raw in lines[data_start - 1:]:
        s = raw.strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            f = int(float(parts[0]))
            v = float(parts[1])
        except ValueError:
            continue
        factor.append(f)
        value.append(v)

    return {
        "factor": factor,
        "value": value,
        "cert_f": cert_f,
        "cert_ssb": cert_ssb,
        "cert_ssw": cert_ssw,
        "data_start_line": data_start,
    }


def compute_fisher(factor, value):
    arr_f = np.asarray(factor)
    arr_v = np.asarray(value, dtype=float)
    groups = [arr_v[arr_f == lvl] for lvl in sorted(set(arr_f))]
    res = stats.f_oneway(*groups)
    f_stat = float(res.statistic)
    p_value = float(res.pvalue)

    # SS via the Welford-style decomposition (matches the TS impl).
    n_total = len(arr_v)
    k = len(groups)
    grand = float(arr_v.mean())
    ssb = float(sum(len(g) * (g.mean() - grand) ** 2 for g in groups))
    ssw = float(sum(((g - g.mean()) ** 2).sum() for g in groups))
    df_bn = k - 1
    df_wn = n_total - k
    df_tot = n_total - 1
    ss_tot = ssb + ssw
    f_critical = float(stats.f.isf(0.05, df_bn, df_wn))
    return {
        "fStat": f_stat,
        "pValue": p_value,
        "ssBn": ssb,
        "ssWn": ssw,
        "ssTot": ss_tot,
        "dfBn": df_bn,
        "dfWn": df_wn,
        "dfTot": df_tot,
        "msBn": ssb / df_bn,
        "msWn": ssw / df_wn,
        "fCritical": f_critical,
    }


def compute_welch(factor, value):
    arr_f = np.asarray(factor)
    arr_v = np.asarray(value, dtype=float)
    groups = [arr_v[arr_f == lvl] for lvl in sorted(set(arr_f))]
    # scipy.stats.f_oneway with equal_var=False is the Welch W-test.
    res = stats.f_oneway(*groups, equal_var=False)
    f_stat = float(res.statistic)
    p_value = float(res.pvalue)

    # Reproduce Welch's df via the formula (matches the TS impl).
    sizes = np.array([len(g) for g in groups], dtype=float)
    variances = np.array([g.var(ddof=1) for g in groups], dtype=float)
    weights = sizes / variances
    W = weights.sum()
    lam = float((((1 - weights / W) ** 2) / (sizes - 1)).sum())
    k = len(groups)
    df_bn = k - 1
    df_wn = (k * k - 1) / (3 * lam)
    f_critical = float(stats.f.isf(0.05, df_bn, df_wn))
    return {
        "fStat": f_stat,
        "pValue": p_value,
        "dfBn": df_bn,
        "dfWn": df_wn,
        "fCritical": f_critical,
    }


# 15-point simple validation dataset (mirrors the existing TS Correctness test).
SIMPLE = {
    "factor": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3],
    "value":  [9, 12, 4, 8, 7, 4, 6, 8, 2, 10, 1, 3, 4, 5, 2],
}


def fmt_num(x: float) -> str:
    if x is None:
        return "null"
    if not np.isfinite(x):
        return "Number.NaN" if np.isnan(x) else ("Number.POSITIVE_INFINITY" if x > 0 else "Number.NEGATIVE_INFINITY")
    # Full-precision Python repr; TS round-trip preserves the same double.
    s = repr(float(x))
    return s


def fmt_int_list(xs):
    return "[" + ", ".join(str(int(v)) for v in xs) + "]"


def fmt_float_list(xs):
    chunks = []
    for i, v in enumerate(xs):
        chunks.append(fmt_num(v))
    return "[" + ", ".join(chunks) + "]"


def emit_fisher(prefix: str, payload: dict, tol: dict) -> str:
    parts = [
        f"  fStat: {fmt_num(payload['fStat'])},",
        f"  pValue: {fmt_num(payload['pValue'])},",
        f"  ssBn: {fmt_num(payload['ssBn'])},",
        f"  ssWn: {fmt_num(payload['ssWn'])},",
        f"  ssTot: {fmt_num(payload['ssTot'])},",
        f"  dfBn: {payload['dfBn']},",
        f"  dfWn: {payload['dfWn']},",
        f"  dfTot: {payload['dfTot']},",
        f"  msBn: {fmt_num(payload['msBn'])},",
        f"  msWn: {fmt_num(payload['msWn'])},",
        f"  fCritical: {fmt_num(payload['fCritical'])},",
        f"  tol: {{rtol: {fmt_num(tol['rtol'])}, atol: {fmt_num(tol['atol'])}}},",
    ]
    return "{\n" + "\n".join(prefix + p for p in parts) + "\n" + prefix + "}"


def emit_welch(prefix: str, payload: dict, tol: dict) -> str:
    parts = [
        f"  fStat: {fmt_num(payload['fStat'])},",
        f"  pValue: {fmt_num(payload['pValue'])},",
        f"  dfBn: {payload['dfBn']},",
        f"  dfWn: {fmt_num(payload['dfWn'])},",
        f"  fCritical: {fmt_num(payload['fCritical'])},",
        f"  tol: {{rtol: {fmt_num(tol['rtol'])}, atol: {fmt_num(tol['atol'])}}},",
    ]
    return "{\n" + "\n".join(prefix + p for p in parts) + "\n" + prefix + "}"


def main():
    out_path = Path(__file__).resolve().parent.parent / "src" / "tests" / "anova-fixtures.ts"

    blocks: list[str] = []

    # validation_simple
    fisher_s = compute_fisher(SIMPLE["factor"], SIMPLE["value"])
    welch_s = compute_welch(SIMPLE["factor"], SIMPLE["value"])
    simple_block = [
        "  validation_simple: {",
        f"    factor: {fmt_int_list(SIMPLE['factor'])},",
        f"    value: {fmt_int_list(SIMPLE['value'])},",
        "    fisher: " + emit_fisher("    ", fisher_s, SIMPLE_TOLERANCE) + ",",
        "    welch: " + emit_welch("    ", welch_s, SIMPLE_TOLERANCE) + ",",
        "  },",
    ]
    blocks.append("\n".join(simple_block))

    # NIST datasets
    nist_blocks = ["  nist: {"]
    for name in DATASETS:
        print(f"fetch {name}...", file=sys.stderr)
        text = fetch(name)
        parsed = parse_nist_dat(text)
        fisher = compute_fisher(parsed["factor"], parsed["value"])
        welch = compute_welch(parsed["factor"], parsed["value"])
        tol_f = TOLERANCE[name]["fisher"]
        tol_w = TOLERANCE[name]["welch"]
        nist_blocks.append(f"    {name}: {{")
        nist_blocks.append(f"      factor: {fmt_int_list(parsed['factor'])},")
        nist_blocks.append(f"      value: {fmt_float_list(parsed['value'])},")
        nist_blocks.append(f"      certifiedF: {fmt_num(parsed['cert_f'])},")
        nist_blocks.append(f"      certifiedSsBn: {fmt_num(parsed['cert_ssb'])},")
        nist_blocks.append(f"      certifiedSsWn: {fmt_num(parsed['cert_ssw'])},")
        nist_blocks.append("      fisher: " + emit_fisher("      ", fisher, tol_f) + ",")
        nist_blocks.append("      welch: " + emit_welch("      ", welch, tol_w) + ",")
        nist_blocks.append("    },")
    nist_blocks.append("  },")
    blocks.append("\n".join(nist_blocks))

    header = """\
/* AUTO-GENERATED by tools/generate-anova-fixtures.py — do not edit by hand.
 *
 * Sources:
 *   - scipy.stats.f_oneway  (Fisher and Welch W-test)
 *   - NIST StRD ANOVA datasets, https://www.itl.nist.gov/div898/strd/anova/
 *
 * To regenerate, run `python tools/generate-anova-fixtures.py` from packages/EDA.
 */

/* eslint-disable max-len */

export type Tol = {rtol: number; atol: number};

export type FisherFixture = {
  fStat: number;
  pValue: number;
  ssBn: number;
  ssWn: number;
  ssTot: number;
  dfBn: number;
  dfWn: number;
  dfTot: number;
  msBn: number;
  msWn: number;
  fCritical: number;
  tol: Tol;
};

export type WelchFixture = {
  fStat: number;
  pValue: number;
  dfBn: number;
  dfWn: number;
  fCritical: number;
  tol: Tol;
};

export type SimpleDataset = {
  factor: number[];
  value: number[];
  fisher: FisherFixture;
  welch: WelchFixture;
};

export type NistDataset = {
  factor: number[];
  value: number[];
  certifiedF: number | null;
  certifiedSsBn: number | null;
  certifiedSsWn: number | null;
  fisher: FisherFixture;
  welch: WelchFixture;
};

export const FIXTURES = {
"""

    footer = "} as const;\n"

    out = header + "\n".join(blocks) + "\n" + footer
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(out, encoding="utf-8")
    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
