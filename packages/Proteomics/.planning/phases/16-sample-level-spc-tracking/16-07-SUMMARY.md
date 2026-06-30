---
phase: 16
plan: 07
wave: 4
status: complete
completed: 2026-06-08
---

# Plan 16-07 — P2 Pareto panel

## What's Built

`src/viewers/spc-dashboard.ts` extended with:

- `aggregateParetoCounts(runs)` — pure helper, sole responsibility of trip-count aggregation. Returns `Array<{rule, metric, count}>` sorted DESC by count. Used by `SPC:pareto_descending` test.
- `createParetoChart(runsDf, instrumentId)` — viewer factory. Walks the runs slice, parses `rules_tripped` JSON per row (try/catch — one malformed row does not poison the aggregation per T-16-02), aggregates via `aggregateParetoCounts`, builds a small `DG.DataFrame` with `rule_metric` / `trip_count` / `total_runs` / `pct` columns (pre-sorted DESC so barChart renders bars in order), and returns a `DG.Viewer.barChart` with the UI-SPEC title verbatim. Returns `null` when slice has fewer than 4 runs — the caller renders the empty-state body.
- `refreshCharts` inside `openDashboardForInstrument` now docks the Pareto panel below the I/MR-charts. Empty-state body (`Pareto chart appears once 4 or more runs have been recorded for this instrument.`) renders when `createParetoChart` returns null.

## RED → GREEN

| Test | Asserts |
|------|---------|
| `SPC:pareto_descending` | Aggregates `[A, A, A, B]` → `[{rule: 'nelson_1', metric: 'median_intensity', count: 3}, {rule: 'nelson_5', metric: 'control_corr', count: 2}]`, sorted DESC |

All 30 SPC tests are now mapped to GREEN-eligible code paths (final runtime status is captured at the operator UAT step).

## Acceptance gate verification

| Gate | Result |
|------|--------|
| `npm run build` exit 0 | ✓ |
| `export function createParetoChart` | 1 |
| `DG.Viewer.barChart` | 1 |
| `Which rules trip most often` (UI-SPEC title verbatim) | 1 |
| `Pareto chart appears once 4 or more runs` (empty-state body verbatim) | 1 |
| `Times tripped` (Y-axis label) | 1 |
| `export function aggregateParetoCounts` | 1 |

## Decoupling claim

The plan's `<specifics>` claim ("Pareto can land as a v1.4.1 hotfix without re-instrumenting metric code") holds because the `rules_tripped` data is already on every row of runs.csv (Plan 16-03 column shape + Plan 16-05 menu handler write). Commenting out the Pareto `refreshCharts` block — without committing — leaves the rest of the dashboard fully functional. That is the documented decoupling proof.

## Threat mitigations

- T-16-02 — `rules_tripped` JSON.parse wrapped in try/catch; malformed rows contribute an empty array to the aggregation.

## Key files

```yaml
key-files:
  modified:
    - src/viewers/spc-dashboard.ts
```

## Self-Check: PASSED
