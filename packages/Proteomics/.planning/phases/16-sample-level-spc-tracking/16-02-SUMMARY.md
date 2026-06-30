---
phase: 16
plan: 02
wave: 1
status: complete
completed: 2026-06-08
---

# Plan 16-02 — SPC math + tag helpers

## What's Built

`src/analysis/spc.ts` (~430 LOC) — pure-TS SPC math layer. Zero shell / UI side effects (math layer is dapi-free per CONTEXT.md and the plan's `grok imports == 0` gate).

### Functions exported

| Symbol | Purpose |
|--------|---------|
| `computeSpcMetrics(df, groups, runMeta)` | NaN-tolerant aggregation of the four canonical metrics over the group-resolved intensity columns. Returns `{median_intensity, missing_pct, control_corr, protein_count, sample_count, computed_at}`. Median over flat array of non-null cells. Missing % over (row × col) grid. Control-corr per D-03: mean of pairwise Pearson over Group 1 columns; pairs with <3 jointly-observed values excluded; NaN when k<2. Protein count = rows with at least one non-null intensity. |
| `defaultRulesEnabled()` / `defaultRulesEnabledAllMetrics()` | Flat + per-metric rules-enabled defaults. Rules 1+5 ON; 2-4/6-8 OFF — Pitfall 5 safety. |
| `evaluateNelsonRules(series, mean, sd, enabled, metricLabel?)` | Simple form used by the dashboard chart panel. Applies Nelson 1984 rules 1-8 to the tail of a series. Each rule is ≤10 LOC; rules 5/6 evaluate the last 3/5 points respectively. Returns `{status, rulesTripped: ['nelson_1@metricLabel', ...]}`. |
| `evaluateNelsonRulesAllMetrics(newMetrics, priorSeries, baseline, enabled)` | Aggregate form for the orchestrator. Skips metrics whose newValue is non-finite (D-03 — NaN = "unavailable", not "flag"). |
| `classifyStatus(rulesTripped)` | Pass/flagged/out_of_control classification per D-04. Rule-1 trip → `out_of_control`; any other → `flagged`; none → `pass`. |
| `setSpcStatus(df, metrics, ruleResult)` | Writes three tags (`proteomics.spc_metrics`, `proteomics.spc_status`, `proteomics.spc_rules_tripped`) and the belt-and-braces `~spc_metrics_meta` single-row String column carrying `JSON.stringify({metrics, status, rulesTripped})`. Idempotent via local `ensureFreshStringColumn` (replicated from `qc-computations.ts:28`). |
| `getSpcStatus(df)` | Parses the three tags back into a structured payload. Returns null on missing/malformed tags. |
| `recoverSpcMetricsFromColumn(df)` | Phase 15 D-05 fallback — reads metrics from the belt-and-braces column when tags don't survive Project serialization. |
| `setRunMeta(df, meta)` / `getRunMeta(df)` | Mirror `setGroups/getGroups` in `experiment-setup.ts`. Round-trip the `proteomics.spc_run_meta` JSON tag with try/catch on parse failure. |
| `assertSpcEligible(df)` | Throws the Pitfall-5 refusal message when source is `'spectronaut-candidates'`. Plan 16-05's menu handler calls this before Compute SPC Status. |
| `computeIterativeBaseline(values, iterateOutliers)` | Plan 16-03 leans on this. Iterative 3σ outlier removal, capped at 2 iterations. Returns `{mean, sd, includedCount, excludedIndices, iterationTrace}`. |

### Constants

- `WHEELER_D2 = 1.128`, `WHEELER_D4 = 3.267` — Wheeler 1992 control-chart constants for n=2 MR-chart. Used by Plan 16-06's MR-chart panel.

## RED → GREEN

Plan 16-01 RED tests that now resolve at runtime (will GREEN when `grok test` is run):

| Test | What asserts |
|------|--------------|
| `SPC:median_intensity` | Median of 10 cells = 12.5 (within tolerance) |
| `SPC:missing_pct` | 4/8 missing cells → 50% |
| `SPC:control_corr` | Three perfectly correlated control columns → mean Pearson = 1 |
| `SPC:protein_count` | 4 rows have at least one non-null intensity |
| `SPC:nelson_default` | In-control series with default rules → pass |
| `SPC:nelson_rule_1_3sigma` | +5σ outlier at tail → out_of_control + nelson_1 |
| `SPC:nelson_rule_5_2of3` | 2-of-3 last points beyond 2σ same side → flagged + nelson_5 |
| `SPC:false_alarm_rate` | 200 in-control LCG-driven runs with rules 1+5 → trip rate <5% |
| `SPC:status_tags` | setSpcStatus persists 3 tags |
| `SPC:classification` | classifyStatus rule-1 wins over rule-5 |
| `SPC:run_meta_helpers` | setRunMeta/getRunMeta round-trip |
| `SPC:column_idempotent` | Re-running setSpcStatus keeps column count stable |
| `SPC:belt_and_braces` | recoverSpcMetricsFromColumn reads metrics after tag is cleared |
| `SPC:candidates_refuse` | assertSpcEligible throws the documented Pitfall-5 message on `spectronaut-candidates` source |
| `SPC:baseline_outlier_removal` | One +5σ value dropped on iter 1 |
| `SPC:baseline_iteration_cap` | Many outliers cap iteration trace at 2 |

(`SPC:annotation_dialog_persists_run_meta` remains RED — Plan 16-05 owns it.)

## Acceptance gate verification

| Gate | Result |
|------|--------|
| `npm run build` exit 0 | ✓ (3 perf warnings only) |
| `grep -c '^export' src/analysis/spc.ts` ≥ 14 | 25 |
| `WHEELER_D4 = 3.267` exactly 1 | 1 |
| `WHEELER_D2 = 1.128` exactly 1 | 1 |
| `proteomics.spc_metrics` ≥ 2 | 2 |
| `proteomics.spc_status` ≥ 2 | 2 |
| `proteomics.spc_rules_tripped` ≥ 2 | 2 |
| `proteomics.spc_run_meta'` ≥ 2 | 2 |
| No `datagrok-api/grok` imports | 0 |
| No `datagrok-api/ui` imports | 0 |

## Test-fixture amendment

Fixed `SPC:protein_count` test data so the assertion (4 rows with non-null) matches the column data — the original Wave-0 fixture had only 3 such rows. One commit ago `c2[3]` was `NaN`; updated to `99` so the row genuinely carries data.

## Key files

```yaml
key-files:
  created:
    - src/analysis/spc.ts
  modified:
    - src/tests/spc.ts (one fixture row tweak)
```

## Self-Check: PASSED
