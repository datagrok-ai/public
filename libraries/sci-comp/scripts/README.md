# Fixture Generators

Generates JSON fixtures and source datasets used by the TypeScript test
harnesses. All scripts are run rarely — the resulting CSV/JSON files are
committed to git, so CI does not need Python, R, or any external repo.

| Script | Output | Tooling |
|---|---|---|
| `generate-fixtures.py` | `src/stats/__tests__/fixtures/*.json` | Python + scipy |
| `normalize-jonckheere-fixtures.py` | `src/stats/__tests__/fixtures/jonckheere-*.json` | Python (no extras) |
| `generate-dataset-03-rat-simple.R` | `src/nca/__tests__/datasets/03_rat_simple.csv` | R |
| `generate-nca-fixtures.R` | `src/nca/__tests__/fixtures/*.json` | R + PKNCA |

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

### Requirements

Plain Python 3.9+ — no extra packages. A local checkout of
`MatthewCorney/jonckheere_terpstra`.

### How to run

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

## NCA Fixtures (R + PKNCA)

Reference values for the NCA Phase 1 validation gate (see
`docs/nca_development_plan_v2.md`, Tasks 1.7.5 and 1.9.5). Two scripts:

- `generate-dataset-03-rat-simple.R` — produces the synthetic rat oral PK
  CSV in `src/nca/__tests__/datasets/03_rat_simple.csv`. Datasets 01
  (Theoph) and 02 (Indometh) are exported one-time from the R `datasets`
  package by hand (see `src/nca/__tests__/datasets/README.md`).
- `generate-nca-fixtures.R` — reads all three CSVs, runs PKNCA with the
  pinned configuration, and writes JSON reference values to
  `src/nca/__tests__/fixtures/*.json`.

### Pinned versions (Phase 0)

| Tool | Version |
|---|---|
| R | 4.6.0 (ucrt build, 2026-04-24) |
| PKNCA | 0.12.1 |
| jsonlite | 2.0.0 |

### PKNCA configuration (mirrored in `NcaRules` defaults on TS side)

| Option | Value | Notes |
|---|---|---|
| `auc.method` | `lin up/log down` | PKNCA 0.12 naming; was `linear up log down` pre-0.10 |
| `min.points` | 3 | Minimum points for lambda_z best-fit |
| `min.r.squared` | 0.85 | Minimum adjusted R² to accept lambda_z |
| `exclude.cmax` | TRUE | Cmax point excluded from lambda_z candidates |
| `min.span.ratio` | 2 | PKNCA default — lambda_z span must cover ≥ 2 × t½ |
| `max.aucinf.pext` | 50 | % AUC extrapolated — hard error threshold |
| `conc.blq` | `{first: keep, middle: drop, last: drop}` | BLQ handling |

### PKNCA 0.12 quirks the script handles

- **Result column renamed**: `PPSTRES` → `PPORRES`.
- **AUC interval requires `start ≥ first measurable time`**. Datasets without
  t=0 (Indometh, rat_simple) need a t=0 point inserted to capture the
  `[0, t_first]` portion of AUC. Insertion strategy depends on route:
  - `extravascular`: insert (0, 0) (pre-dose conc is 0 by convention).
  - `iv-bolus`: two-pass — pass 1 estimates c0 via PKNCA's `c0` parameter
    (NA at t=0); pass 2 inserts (0, c0) and runs full NCA. Cmax/Tmax are
    then post-processed to be the OBSERVED peak (not the inserted c0).
- **`NULL` inside R lists serialises as `{}` in jsonlite**. Use `unbox(NA)`
  combined with `toJSON(na = "null")` to get JSON `null`.
- **`min.span.ratio = 0` rejected**. Use `2` (PKNCA default) or higher.

### Usage

Set personal library env var first (one-time):

```powershell
[Environment]::SetEnvironmentVariable('R_LIBS_USER', "$env:LOCALAPPDATA\R\win-library\4.6", 'User')
```

Then from `libraries/sci-comp/`:

```powershell
# Regenerate dataset 03 (rare — only if you want different noise/seed):
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" scripts/generate-dataset-03-rat-simple.R

# Regenerate all NCA fixtures from the three CSVs:
& "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" scripts/generate-nca-fixtures.R
```

Output: three JSON files under `src/nca/__tests__/fixtures/`:
- `01_theoph.json` — 12 profiles
- `02_indometh.json` — 6 profiles
- `03_rat_simple.json` — 8 profiles

### JSON shape

Per dataset:

```json
{
  "dataset": "<name>",
  "pknca_version": "0.12.1",
  "config": { "auc_method": "...", "min_points": 3, ... },
  "dataset_meta": { "dose": { "value": ..., "unit": "..." }, "route": "..." },
  "profiles": [
    {
      "profile_key": { "subject": "...", "route": "..." },
      "parameters": {
        "cmax": ..., "tmax": ..., "auclast": ..., "aucinf": ...,
        "pct_aucextrap": ..., "lambda_z": ..., "half_life": ...,
        "cl": ..., "vz": ...
      },
      "provenance": {
        "lambda_z_n_points": ..., "lambda_z_r_squared": ...,
        "lambda_z_adj_r_squared": ..., "lambda_z_time_first": ...,
        "lambda_z_time_last": ..., "clast_obs": ...,
        "c0_extrapolated": ...
      }
    }
  ]
}
```

`cl` and `vz` are PKNCA-style universal keys. Their physical interpretation
(`CL` vs `CL/F`, `Vz` vs `Vz/F`) follows from the `route` field — see
`src/nca/__tests__/datasets/README.md` and Variant 1 convention in
`docs/nca_development_plan_v2.md` §6.
