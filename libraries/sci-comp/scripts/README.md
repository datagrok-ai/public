# Fixture Generator

Generates JSON fixtures for the TypeScript stats port from the original Python
validation suite. Run once (or on scipy version bumps) — the resulting JSON
files in `src/stats/__tests__/fixtures/` are committed to git, so CI does not
need Python or scipy.

## Prerequisites

Conda environment `SEND-TEST` with `numpy` and `scipy>=1.11`.

## Usage

```bash
# From repo root or anywhere — the script uses absolute paths internally.
"$USERPROFILE/anaconda3/envs/SEND-TEST/python.exe" \
  c:/Users/vmaka/Datagrok/public-copy/public/libraries/sci-comp/scripts/generate-fixtures.py
```

Output: `src/stats/__tests__/fixtures/*.json`. One file per Python
validation script.

## Provenance

Each generated JSON includes a `metadata` block with:
- `source`: name of the original `validate_*.py`
- `scipy_version`, `numpy_version`, `python_version`
- `tolerances`: default tolerances inherited from `validate_helpers.py`

When tolerances are overridden per-case in the original Python code, the
override is recorded in the case-level `tolerance` field.
