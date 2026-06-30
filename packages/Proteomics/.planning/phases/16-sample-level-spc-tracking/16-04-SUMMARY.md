---
phase: 16
plan: 04
wave: 1
status: complete
completed: 2026-06-08
---

# Plan 16-04 — Spectronaut PG-report run-meta seed

## What's Built

`src/parsers/spectronaut-parser.ts` augmented in lockstep on both paths:

### Text path (`pivotSpectronaut`)

- Extended `PivotResult` interface with two optional fields: `earliestRunDate?` (ISO-8601 or null) and `firstInstrumentMethod?` (verbatim cell or null).
- Single-pass extraction inside the existing iteration over `longDf` — Plan Perf Decision Option A. No second file pass; satisfies the multi-GB Spectronaut PG-report constraint (memory `project_proteomics_large_spectronaut_aggregation`).
- Date parsing via the native `Date.parse` (vs adding `dayjs` — no existing usage in this package); invalid dates skip rather than poison the earliest-observed comparison.

### Streaming path (`parseSpectronautStream`)

- Closure-scoped accumulators `earliestRunDateMs` / `earliestRunDateIso` / `firstInstrumentMethod`.
- Optional `runDateI` + `instrMethodI` index lookups inside `parseHeader`; both default to `-1` when the columns are absent.
- `expectedFields` now includes the two new indices in its `Math.max(...)` call (only when `>= 0` because `-1` is harmless for `Math.max(... , -1)` + 1 — already 0-floored by the protI/condI minimum guarantee).
- `handleFields` extracts the date and instrument method inline on every kept row.
- Aggregator handoff: `aggToPivotResult(agg)` returns a `PivotResult`; the streaming path overwrites `earliestRunDate` / `firstInstrumentMethod` on that result with the closure values before invoking `finalizeSpectronaut`.

### Shared helper

- New private `applySpcRunMetaSeed(df, result)` — single point that writes `proteomics.spc_run_meta_seed` JSON tag (`{instrument_id, acquisition_datetime}`). Called from both the streaming and text exit paths. The tag is NOT set when both fields are absent — backward compatibility for v1.3 inputs.

## RED → GREEN

| Test | Asserts |
|------|---------|
| `SPC:spectronaut_seed` | Feeding a synthetic header with both columns yields `proteomics.spc_run_meta_seed = {"instrument_id": "QExactive-01", "acquisition_datetime": "2026-07-10T00:00:00.000Z"}` |

## Backward compatibility

`REQUIRED_COLUMNS` is unchanged (`['R.Condition', 'R.Replicate', 'PG.ProteinGroups']`). A v1.3 PG export without `R.RunDate` / `R.InstrumentMethod` produces a df with NO seed tag — Plan 16-05's dialog will fall back to empty inputs. Existing Spectronaut tests under `src/tests/spectronaut-parser.ts` / `src/tests/parsers.ts` are not modified and remain GREEN.

## Acceptance gate verification

| Gate | Result |
|------|--------|
| `npm run build` exit 0 | ✓ |
| `R.RunDate` references | 5 |
| `R.InstrumentMethod` references | 4 |
| `proteomics.spc_run_meta_seed` references | 3 |
| `REQUIRED_COLUMNS` does NOT include new columns | confirmed (3 entries unchanged) |

## Key files

```yaml
key-files:
  modified:
    - src/parsers/spectronaut-parser.ts
```

## Self-Check: PASSED
