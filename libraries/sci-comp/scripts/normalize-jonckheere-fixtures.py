"""Normalise Jonckheere-Terpstra reference fixtures from the external Python
repo (https://github.com/MatthewCorney/jonckheere_terpstra) into the format
expected by `__tests__/helpers.ts::loadFixture`.

Each source JSON is a flat list of records with `x`, `g`, `alt`, `nperm`,
`continuity`, `statistic` (J), `zstat`, `p_value`, `significant`. We emit
JSON files at `src/stats/__tests__/fixtures/jonckheere-{slug}.json` with the
shape `{metadata, cases}`.

This is a one-shot normaliser — fixtures are committed to git, so it is
rarely re-run (only after a refresh of the upstream reference data).

Usage:
    # Default: looks for the external repo at the path below, but any
    # checkout of MatthewCorney/jonckheere_terpstra can be supplied via
    # the JT_REPO env var or the --repo CLI argument.
    python scripts/normalize-jonckheere-fixtures.py
    JT_REPO=/path/to/jonckheere_terpstra python scripts/normalize-jonckheere-fixtures.py
    python scripts/normalize-jonckheere-fixtures.py --repo /path/to/jonckheere_terpstra
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

HERE = Path(__file__).resolve().parent
LIB_ROOT = HERE.parent
FIXTURES_DIR = LIB_ROOT / "src" / "stats" / "__tests__" / "fixtures"

# Developer-machine default; override with JT_REPO env var or --repo flag.
DEFAULT_REPO = Path(r"C:/Users/vmaka/Datagrok/NumericLibsPy/jonckheere_terpstra")


def _norm_alt(alt: str) -> str:
    """Map the Python repo's alternative codes to our `Alternative` type."""
    if alt == "two_sided":
        return "two-sided"
    if alt in ("increasing", "decreasing"):
        return alt
    raise ValueError(f"Unknown alternative: {alt!r}")


def _normalise(
    repo_root: Path,
    src: Path,
    method: str,
    tolerances: dict[str, float],
    name_prefix: str,
    *,
    continuity_override: bool | None = None,
) -> dict[str, Any]:
    """Read the external JSON, wrap each record into our fixture-case shape."""
    raw: list[dict[str, Any]] = json.loads(src.read_text(encoding="utf-8"))
    cases: list[dict[str, Any]] = []
    for i, rec in enumerate(raw):
        nperm = rec.get("nperm")
        continuity = continuity_override if continuity_override is not None else (
            bool(rec.get("continuity")) if isinstance(rec.get("continuity"), bool) else True
        )
        cases.append({
            "name": f"{name_prefix}#{i}: alt={rec['alt']}, "
                    f"k={rec.get('group_number', '?')}, n={rec.get('group_size', '?')}, "
                    f"id={rec.get('id', '?')}"
                    + (f", nperm={nperm}" if nperm else ""),
            "category": method,
            "inputs": {
                "x": [float(v) for v in rec["x"]],
                "g": [int(v) for v in rec["g"]],
                "alternative": _norm_alt(rec["alt"]),
                "method": method,
                "continuity": continuity,
                "nperm": int(nperm) if isinstance(nperm, int) else None,
            },
            "expected": {
                "j_statistic": float(rec["statistic"]),
                "z_statistic": (float(rec["zstat"])
                                if isinstance(rec.get("zstat"), (int, float))
                                else None),
                "p_value": float(rec["p_value"]),
                "significant": bool(rec.get("significant", False)),
            },
        })
    return {
        "metadata": {
            "source": str(src.relative_to(repo_root)),
            "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "method": method,
            "tolerances": tolerances,
        },
        "cases": cases,
    }


def write(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {path.relative_to(LIB_ROOT)} ({len(payload['cases'])} cases)")


def resolve_repo() -> Path:
    """Resolve the external repo root from CLI / env / default, in that order."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo",
        type=Path,
        default=None,
        help="Path to a checkout of MatthewCorney/jonckheere_terpstra. "
             "Falls back to $JT_REPO, then to a hardcoded developer default.",
    )
    args = parser.parse_args()
    return args.repo or Path(os.environ.get("JT_REPO", str(DEFAULT_REPO)))


def main() -> int:
    repo_root = resolve_repo()
    tests_root = repo_root / "tests"
    if not tests_root.exists():
        print(
            f"external fixture root not found: {tests_root}\n"
            f"  Pass --repo /path/to/jonckheere_terpstra or set $JT_REPO.",
            file=sys.stderr,
        )
        return 1
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)

    # NOTE on tolerances: `expectClose` uses `diff >= tol`, so tol=0 always
    # fails. We use the midrank convention for J (matches clinfun / PMCMR /
    # DescTools / SAS PROC FREQ), so untied data must match R exactly →
    # use 0.5 (integer match window). R fixtures show occasional 0.5–1.0
    # offsets on tied cases (R serialises J to integer, dropping the .5
    # half-step) → use 1.5 there. z and p still match tightly.

    # 1. regressionpack approximate. Tolerance matches the upstream pytest
    # suite (rel=0.10, abs=0.02 → use abs ~0.05 to cover relative noise on
    # mid-range p-values). regressionpack uses a slightly different variance
    # formula than Lehmann §6.2, hence the slack.
    write(
        FIXTURES_DIR / "jonckheere-rp-approximate.json",
        _normalise(
            repo_root,
            tests_root / "python_generation" / "jonckheere_test_results_regression_pack.json",
            method="approximate",
            tolerances={"p_value": 0.05, "z_statistic": 0.5, "j_statistic": 0.5},
            name_prefix="rp-approx",
            continuity_override=True,  # regressionpack always applies continuity
        ),
    )

    # 2. regressionpack exact (no ties in fixture data → J must match exactly)
    write(
        FIXTURES_DIR / "jonckheere-rp-exact.json",
        _normalise(
            repo_root,
            tests_root / "python_generation" / "exact_jonckheere_test_results_regression_pack.json",
            method="exact",
            tolerances={"p_value": 1e-3, "j_statistic": 0.5},
            name_prefix="rp-exact",
        ),
    )

    # 3. R clinfun approximate (cross-implementation). j tolerance 1.5
    # covers R's integer-rounded J on tied cases.
    write(
        FIXTURES_DIR / "jonckheere-clinfun.json",
        _normalise(
            repo_root,
            tests_root / "R_generation" / "jonckheere_test_results_clinfun.json",
            method="approximate",
            tolerances={"p_value": 0.05, "j_statistic": 1.5},
            name_prefix="clinfun",
            continuity_override=True,
        ),
    )

    # 4. R PMCMRplus permutation (cross-implementation; loose tolerance —
    # PRNG differs, so this is a Monte Carlo vs Monte Carlo comparison).
    # Upstream pytest uses rel=0.15, abs=0.02 → use abs ~0.10 so that
    # noisy mid-range p-values pass. j tolerance 1.5 same rationale as #3.
    write(
        FIXTURES_DIR / "jonckheere-pmcmr.json",
        _normalise(
            repo_root,
            tests_root / "R_generation" / "jonckheere_test_results_PMCMRplus.json",
            method="permutation",
            tolerances={"p_value": 0.10, "j_statistic": 1.5},
            name_prefix="pmcmr",
        ),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
