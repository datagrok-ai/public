---
phase: 16
plan: 01
wave: 0
status: complete
completed: 2026-06-08
---

# Plan 16-01 — Wave 0 RED scaffold + formula-line spike

## What's Built

- `src/tests/spc-formula-lines-spike.ts` — automated spike that calls `viewer.meta.formulaLines.addLine({style:'dashed', ...})` on a `DG.Viewer.lineChart` over a 5-row datetime/float DataFrame, writes the verbatim observation to `System:AppData/Proteomics/spike/spike-result.md`, and ALWAYS green-passes (`expect(true, true)`) — the spike's purpose is to ANSWER the question, not to gate. Two test() entries: `SPC:spike_formula_line_dashed_on_datetime` (single line) + `SPC:spike_three_lines_ucl_cl_lcl` (full UCL/CL/LCL stack).
- `src/tests/spc.ts` — 30 RED test stubs under `category('SPC', ...)` mapped one-to-one to VALIDATION.md rows and RESEARCH.md "Phase Requirements → Test Map". Test ownership split (run-meta helper round-trip vs. annotation-dialog persistence) is implemented as two separate tests (`SPC:run_meta_helpers` for Plan 16-02; `SPC:annotation_dialog_persists_run_meta` for Plan 16-05).
- `src/package-test.ts` — registers both new test files after the existing `publish-roundtrip` import and before `export const _package`.

## Spike outcome

The spike test is automated and self-recording. Wave 0 acceptance only requires that the spike RUN to completion regardless of `PASS`/`FAIL`. The human-verified outcome will be captured in `16-01-SPIKE-RESULT.md` during the checkpoint by running the test on the operator's Datagrok instance and copying `System:AppData/Proteomics/spike/spike-result.md` into the planning folder. Plan 06's `createSpcChartPanel` branches deterministically on this marker:
- `OUTCOME: PASS` → use `style: 'dashed'` on UCL/LCL formula lines.
- `OUTCOME: FAIL` → fall back to the constant-column overlay path (RESEARCH.md Pattern 4 fallback).

## RED scaffold inventory

30 tests under `category('SPC')`:

| SPC ID | Tests | Owner plan |
|--------|-------|------------|
| SPC-01 | `median_intensity`, `missing_pct`, `control_corr`, `protein_count` | 16-02 |
| SPC-03 | `nelson_default`, `nelson_rule_1_3sigma`, `nelson_rule_5_2of3`, `false_alarm_rate` | 16-02 |
| SPC-04 | `status_tags`, `classification` | 16-02 |
| SPC-05 | `baseline_outlier_removal`, `baseline_iteration_cap`, `baseline_roundtrip`, `baseline_rebuild_overwrites`, `rule_toggle_per_instrument` | 16-02 (math), 16-03 (storage) |
| SPC-06 | `run_meta_helpers`, `annotation_dialog_persists_run_meta`, `spectronaut_seed`, `backfill_ordering`, `candidates_refuse` | 16-02, 16-04, 16-05 |
| D-02/Pitfall 8 | `idempotent_upsert`, `storage_bounded`, `schema_version`, `column_idempotent`, `belt_and_braces` | 16-03 (storage), 16-02 (column) |
| SPC-07/08 | `drilldown_resolved`, `drilldown_missing`, `dashboard_renders`, `formula_lines`, `pareto_descending` | 16-06, 16-07 |

Synthetic-fixture helpers (`makeLcg`, `gaussian`, `makeInControlSeries`, `makeOutOfControlSeries`, `makeSyntheticRunsDf`, `makeSyntheticSpectronautHeader`, `fixturePath`, `cleanupFixture`) live as non-exported top-level functions inside `src/tests/spc.ts` per the plan contract.

## Acceptance gate verification

| Gate | Command | Result |
|------|---------|--------|
| Test count ≥ 29 | `grep -c "test(" src/tests/spc.ts` | 30 |
| `test('SPC:` count ≥ 29 | `grep -c "test('SPC:" src/tests/spc.ts` | 30 |
| `SPC:run_meta_helpers` exactly 1 | grep | 1 |
| `SPC:annotation_dialog_persists_run_meta` exactly 1 | grep | 1 |
| Sandbox path coverage ≥ 3 | `grep -c "System:AppData/Proteomics/spc-test/"` | 4 |
| Dashed-style markers ≥ 2 in spike | `grep -c "style.*dashed" src/tests/spc-formula-lines-spike.ts` | 5 |
| `npm run build` exit 0 | `npm run build` | 0 (3 perf warnings only) |
| Test-file imports in package-test.ts | grep | both present |

## Forward-reference convention

Tests that target not-yet-implemented modules (`../analysis/spc`, `../analysis/spc-storage`, `../viewers/spc-dashboard`) use dynamic `await import(/* webpackIgnore: true */ '<path>')` calls with a `// @ts-ignore` annotation. The `webpackIgnore: true` magic comment tells webpack to leave the import as a native browser `import()` — at Wave 0 the modules don't exist so the runtime resolution fails with a deterministic "Failed to fetch dynamically imported module" RED error. Plans 02/03/06 create the destination files; the runtime imports begin resolving as those plans land.

## Key files

```yaml
key-files:
  created:
    - src/tests/spc-formula-lines-spike.ts
    - src/tests/spc.ts
  modified:
    - src/package-test.ts
```

## Wave 0 status

Locked: API contract for `DG.Viewer.lineChart` formula lines on datetime X will be answered by the spike at runtime; 30 RED tests are queued for Waves 1–4 to turn GREEN.

## Self-Check: PASSED
