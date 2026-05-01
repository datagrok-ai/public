# Fixture Generators

Generates JSON fixtures in `src/stats/__tests__/fixtures/` for the TypeScript
stats port. Both scripts are run rarely — the resulting JSON files are
committed to git, so CI does not need Python, scipy, or any external repo.

## `generate-fixtures.py` — SEND-TEST validation suite

Computes reference values live from scipy / numpy. Re-run on scipy version
bumps or when a new test case is added.

### Prerequisites

Conda environment `SEND-TEST` with `numpy` and `scipy>=1.11`.

### Usage

```bash
# From repo root or anywhere — the script uses absolute paths internally.
"$USERPROFILE/anaconda3/envs/SEND-TEST/python.exe" \
  c:/Users/vmaka/Datagrok/public-copy/public/libraries/sci-comp/scripts/generate-fixtures.py
```

Output: `src/stats/__tests__/fixtures/*.json` — one file per Python
validation script.

### Provenance

Each generated JSON includes a `metadata` block with:
- `source`: name of the original `validate_*.py`
- `scipy_version`, `numpy_version`, `python_version`
- `tolerances`: default tolerances inherited from `validate_helpers.py`

When tolerances are overridden per-case in the original Python code, the
override is recorded in the case-level `tolerance` field.

## `normalize-jonckheere-fixtures.py` — Jonckheere-Terpstra reference data

One-shot normaliser that re-shapes the reference fixtures from the external
[MatthewCorney/jonckheere_terpstra](https://github.com/MatthewCorney/jonckheere_terpstra)
repo into our `{metadata, cases}` format. Re-run only when the upstream
reference data is refreshed.

### Prerequisites

Plain Python 3.9+ — no extra packages. A local checkout of
`MatthewCorney/jonckheere_terpstra`.

### Usage

```bash
# Specify the repo path via --repo or the JT_REPO env var. With no args,
# falls back to a developer-machine default that probably won't exist on
# your machine.
python scripts/normalize-jonckheere-fixtures.py --repo /path/to/jonckheere_terpstra
JT_REPO=/path/to/jonckheere_terpstra python scripts/normalize-jonckheere-fixtures.py
```

Output: four files under `src/stats/__tests__/fixtures/`:
- `jonckheere-rp-approximate.json` — Python `regressionpack`, approximate (24 cases).
- `jonckheere-rp-exact.json` — Python `regressionpack`, exact (5 cases).
- `jonckheere-clinfun.json` — R `clinfun::jonckheere.test`, approximate (144 cases).
- `jonckheere-pmcmr.json` — R `PMCMRplus`, permutation (144 cases).

### Tolerances (set in the script, embedded in `metadata.tolerances`)

| Fixture | p_value | j_statistic | Notes |
|---|---|---|---|
| rp-approximate | 0.05 | 0.5 | regressionpack uses a different variance formula → slack mirrors upstream pytest (rel=0.10, abs=0.02). |
| rp-exact | 0.001 | 0.5 | Same exact algorithm — tight match. |
| clinfun | 0.05 | 1.5 | Same midrank J convention; R serialises J to integer, so tied cases drift by ½–1. |
| pmcmr | 0.10 | 1.5 | Monte Carlo vs Monte Carlo with different PRNGs — loose. Same J slack as clinfun. |

`j_statistic` tolerance starts at `0.5` rather than `0` because `expectClose`
in `__tests__/helpers.ts` uses `diff >= tol`, so `tol = 0` always fails on
exact integer matches.
