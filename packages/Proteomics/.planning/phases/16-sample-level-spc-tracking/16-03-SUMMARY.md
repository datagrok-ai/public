---
phase: 16
plan: 03
wave: 1
status: complete
completed: 2026-06-08
---

# Plan 16-03 — SPC AppData storage

## What's Built

`src/analysis/spc-storage.ts` (~400 LOC) — sole owner of the Phase 16 AppData paths. Math layer (`spc.ts`) stays I/O-free; this is the only file in the package that touches `System:AppData/Proteomics/spc/`.

### Public API

| Export | Role |
|--------|------|
| `RUNS_DIR_DEFAULT`, `RUNS_CSV_PATH`, `RUNS_META_PATH`, `RUNS_CSV_SCHEMA_VERSION`, `RUNS_CSV_COLUMNS` | Canonical path + 13-column shape. The columns array IS the Phase 17 additivity contract — Phase 17 appends `campaign_id` to this list and back-fills NULL into existing rows (documented in the file's top JSDoc). |
| `slugifyInstrumentId(raw)` | T-16-01 mitigation. Charset `[A-Za-z0-9._-]+`, 64-char cap, empty → `'unnamed'`. Case preserved (instrument names like `QExactive-01` round-trip). |
| `upsertRun(input, opts?)` / `appendRun` alias | Primary write path. Composite-key upsert on `(instrument_id, acquisition_datetime)`. Second call OVERWRITES the prior row preserving its `run_id`. First-ever invocation pins the schema version into `runs-meta.json`. |
| `loadRuns(instrumentId?, opts?)` | Always sorted by `acquisition_datetime` ASC (Pitfall 7 lever — never insertion order, never `computed_at`). Empty array if `runs.csv` does not exist. |
| `loadBaseline(instrumentId, opts?)` / `saveBaseline(instrumentId, baseline, opts?)` | Per-instrument JSON round-trip. saveBaseline OVERWRITES (no version history — D-04 defers undo to v1.5+). loadBaseline returns null on missing file or JSON parse failure (with `console.warn` — storage layer is shell-free). |
| `readSchemaVersion(opts?)` / `writeSchemaVersion(v, opts?)` / `readRunsMeta(opts?)` / `writeRunsMeta(meta, opts?)` | `runs-meta.json` reader/writer. Idempotent merge — repeated `writeSchemaVersion` calls leave a single file. |
| `iterativeOutlierRemoval(series, maxIter=2)` | Pitfall 6 lever for the storage-side baseline flow (Plan 16-06's modal). Operates on `{run_id, value}` records and returns retained/excluded run_ids + iteration trace. Pathological "every iteration drops one point" inputs hit the cap at exactly 2 iterations. |
| `computeBaselineStats(retainedRows)` | Mean + sample-SD per metric over a retained subset. n-1 denominator for Shewhart 3σ consistency. NaN-tolerant: skips NaN values per metric; returns `{mean: NaN, sd: NaN}` if fewer than 2 finite values. |
| `runsToDataFrame(rows)` | Small adapter used by Plan 16-06's chart panel — turns the parsed CSV rows into a `DG.DataFrame` with a `DATE_TIME` `acquisition_datetime` column. |

### Storage opts

Every public function accepts an optional `{root}` parameter so RED tests can sandbox under `System:AppData/Proteomics/spc-test/<uuid>` (T-16-01-test mitigation) without touching the production `System:AppData/Proteomics/spc/runs.csv`. Production callers omit `root` and the default `RUNS_DIR_DEFAULT` is used.

### CSV encoding choice

The file uses an explicit small CSV writer (single-pass `csvEscape` + header + ordered fields) instead of `DG.DataFrame.toCsv` because the round-trip of mixed `string`/`number`/`null` columns through `DG.DataFrame.fromObjects` muddies type inference on read. The schema is tiny (13 columns, ~104 rows/year) so the simpler hand-written path wins.

## RED → GREEN

| Test | Asserts |
|------|---------|
| `SPC:idempotent_upsert` | Two upserts with same key → 1 row; second-call label wins |
| `SPC:storage_bounded` | 104 weekly cycle inserts collapse to ≤52 per instrument |
| `SPC:schema_version` | `readRunsMeta` returns `{'spc.runs.schema_version': '1'}` after first upsert |
| `SPC:backfill_ordering` | Out-of-order inserts return in ascending `acquisition_datetime` order |
| `SPC:baseline_roundtrip` | saveBaseline → loadBaseline preserves nested structure |
| `SPC:baseline_rebuild_overwrites` | Second saveBaseline wins; only the latest file persists |
| `SPC:rule_toggle_per_instrument` | Custom rules_enabled flags round-trip per instrument |

## Acceptance gate verification

| Gate | Result |
|------|--------|
| `npm run build` exit 0 | ✓ (3 perf warnings only) |
| `RUNS_CSV_PATH` constant | 1 |
| `RUNS_META_PATH` constant | 1 |
| `RUNS_CSV_SCHEMA_VERSION = '1'` | 1 |
| `campaign_id` Phase 17 contract mention | 2 |
| `export function slugifyInstrumentId` | 1 |
| `loadBaseline`/`saveBaseline`/`iterativeOutlierRemoval`/`computeBaselineStats` exports | 1 each |
| `grok.shell` calls | 0 (storage layer is shell-free) |

## Threat mitigations

- **T-16-01** — `slugifyInstrumentId` restricts the instrument-id-to-filename mapping to `[A-Za-z0-9._-]+` and rejects path-traversal characters via the same charset collapse. Empty → `'unnamed'`.
- **T-16-02** — All JSON reads wrapped in try/catch; `runs.csv` reads guarded by `exists()`; numeric coercion via `Number()` is defensive.
- **T-16-03** — Composite-key upsert preserves prior `run_id` on overwrite; Plan 16-05's menu handler emits the "previous status" toast at the user-visible layer.
- **T-16-04** — Single-writer assumption documented in the file's top JSDoc.

## Key files

```yaml
key-files:
  created:
    - src/analysis/spc-storage.ts
```

## Self-Check: PASSED
