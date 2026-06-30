---
phase: 16
plan: 05
wave: 2
status: complete
completed: 2026-06-08
---

# Plan 16-05 — Wire Annotate dialog + top-menu

## What's Built

### `src/analysis/experiment-setup.ts`

- `AnnotationSeed` extended with `instrumentId?` / `acquisitionDatetime?`.
- `seedAnnotationDialogInputs` now resolves the run-meta seed in priority order: persisted `proteomics.spc_run_meta` first, parser-side `proteomics.spc_run_meta_seed` fallback, undefined otherwise.
- `showAnnotationDialog` appends two new inputs AFTER `Group 2`: `Instrument` (`ui.input.string`) and `Acquisition datetime` (`ui.input.date`, date-only precision per RESEARCH gotcha #1 — `ui.input.dateTime` doesn't exist in this js-api). Tooltip copy is VERBATIM from UI-SPEC.
- `applyAnnotation(df, payload)` exported as the programmatic OK-path. Both the live dialog and Plan 16-01's `SPC:annotation_dialog_persists_run_meta` test call this — single source of truth for setGroups + setRunMeta in lockstep. Empty run-meta values are tolerated at annotation time per CONTEXT.md D-01.

### `src/package.ts`

Two new menu items added in the existing decorator pattern (matching `differentialExpression` and `showQcDashboard` shape).

#### `Proteomics | Analyze | Compute SPC Status` (immediate, no `...`)

Gate sequence (all warning copy VERBATIM from UI-SPEC):

1. df present — `Open an analyzed file first, then run Compute SPC Status.`
2. `proteomics.source` ≠ `'spectronaut-candidates'` — `SPC requires per-sample intensities. Re-import this analysis from the Spectronaut PG report (not the Candidates report) to compute SPC.`
3. `getRunMeta(df)` non-empty — `Open Annotate Experiment to set instrument + acquisition datetime.`
4. `getGroups(df)` non-null — `Annotate experimental groups first (Proteomics → Annotate Experiment) — SPC needs Group 1 to compute control-replicate correlation.`
5. Group 1 < 2 cols (info, proceeds) — `Group 1 has fewer than 2 samples — control-replicate correlation can't be computed. The other three metrics will still be recorded.`

Pipeline (`TaskBarProgressIndicator` guarded with try/finally):

1. `computeSpcMetrics(df, groups, runMeta)` → metrics
2. `loadBaseline(runMeta.instrument_id)` → baseline (null on first-ever run)
3. `loadRuns(runMeta.instrument_id)` → priorRuns (used for both prior-status lookup and Nelson rule series)
4. If baseline null: surface `No baseline locked for <instrument_id> yet — recorded as pass. Define a baseline once at least 4 runs have been computed.` and skip rules eval
5. Else: `evaluateNelsonRulesAllMetrics(metrics, priorSeries, baseline.metrics, baseline.rules_enabled)`
6. `setSpcStatus(df, metrics, ruleResult)` writes tags + belt-and-braces column
7. `upsertRun(row)` to runs.csv (idempotent on composite key)
8. Toast selection per status × prior-state:
   - prior row existed → `Updated SPC for <name> (previous status: <prior>).`
   - pass → `SPC computed: <name> — status: pass.`
   - flagged → `SPC computed: <name> — flagged on N rule(s): <list>.`
   - out_of_control → `SPC computed: <name> — OUT OF CONTROL on N rule(s): <list>.`

#### `Proteomics | Visualize | SPC Dashboard...`

Stub. Plan 16-06 replaces the body. Surfaces an info toast pointing at `System:AppData/Proteomics/spc/runs.csv` so the menu is still useful pre-dashboard.

## RED → GREEN

| Test | Asserts |
|------|---------|
| `SPC:annotation_dialog_persists_run_meta` | `applyAnnotation` exists and writes the expected `proteomics.spc_run_meta` JSON |

`SPC:candidates_refuse` was already GREEN (Plan 16-02 owns `assertSpcEligible`). The integrated path in `computeSpcStatus` exercises the same refusal copy from the menu surface.

## Acceptance gate verification

| Gate | Result |
|------|--------|
| `npm run build` exit 0 | ✓ |
| `Proteomics | Analyze | Compute SPC Status` | 1 |
| `Proteomics | Visualize | SPC Dashboard...` | 1 |
| `Open an analyzed file first` | 1 |
| `SPC requires per-sample intensities` | 1 |
| `Open Annotate Experiment to set instrument` | 1 |
| `Annotate experimental groups first` | 3 (existing + new + UI-SPEC body) |
| `Updated SPC for` | 1 |
| `OUT OF CONTROL` | 1 |
| `ui.input.string('Instrument'` | 1 |
| `ui.input.date('Acquisition datetime'` | 1 |
| `ui.input.dateTime` | 0 (RESEARCH gotcha #1) |
| `Free-text identifier for the instrument` tooltip verbatim | 1 |
| `Used to order runs on the SPC trend chart` tooltip verbatim | 1 |

## Threat mitigations

- T-16-01 — `instrument_id` is passed through unchanged to Plan 16-03's `slugifyInstrumentId` for any AppData filename use.
- T-16-03 — `Updated SPC for <name> (previous status: <prior>).` toast makes the idempotent overwrite explicit to the operator.

## Key files

```yaml
key-files:
  modified:
    - src/analysis/experiment-setup.ts
    - src/package.ts
```

## Self-Check: PASSED
