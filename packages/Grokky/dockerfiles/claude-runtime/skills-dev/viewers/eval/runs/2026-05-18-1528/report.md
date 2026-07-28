# Eval report — viewers

Started: 2026-05-18T12:28:58.175Z
Duration: 48.0s
Model: `sonnet`
System prompt size: 26,770 chars
Inlined skills: `datagrok-exec`, `viewers`

## Summary

- **Total prompts:** 20
- **Passed:** 20
- **Failed:** 0
- **Pass rate:** 100.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 8 | 8 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 48.0s wall
- Mean: 8.4s · Median: 7.6s · p95: 18.1s
- Slowest 3: `wrong-tool-trap-chem-space` (18.1s), `add-histogram-split-by-class` (15.7s), `wrong-tool-trap-fromtype-direct` (11.6s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`viewers`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| add-scatter-height-weight | atomic-do | pass | pass | pass | pass | **PASS** | 6.0s |
| add-histogram-age | atomic-do | pass | pass | pass | pass | **PASS** | 6.0s |
| add-line-chart-date-total | atomic-do | pass | pass | pass | pass | **PASS** | 7.4s |
| add-histogram-split-by-class | atomic-do | pass | pass | pass | pass | **PASS** | 15.7s |
| add-bar-chart-count-by-category | atomic-do | pass | pass | pass | pass | **PASS** | 5.5s |
| add-box-plot-activity-by-class | atomic-do | pass | pass | pass | pass | **PASS** | 5.7s |
| add-correlation-plot-numeric | atomic-do | pass | pass | pass | pass | **PASS** | 8.2s |
| close-scatter-plot | atomic-undo | pass | pass | pass | pass | **PASS** | 11.1s |
| close-all-except-grid | atomic-undo | pass | pass | pass | pass | **PASS** | 6.0s |
| remove-histogram | atomic-undo | pass | pass | pass | pass | **PASS** | 7.1s |
| close-every-chart | atomic-undo | pass | pass | pass | pass | **PASS** | 7.5s |
| scatter-mw-logp-colored-regression | combine | pass | pass | pass | pass | **PASS** | 5.6s |
| add-then-close-sequenced | combine | pass | pass | pass | pass | **PASS** | 5.5s |
| replace-scatter-with-heatmap | combine | pass | pass | pass | pass | **PASS** | 8.4s |
| wrong-tool-trap-chem-space | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 18.1s |
| wrong-tool-trap-deprecated-shorthand | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 7.6s |
| wrong-tool-trap-fromtype-direct | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 11.6s |
| reset-to-just-grid | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 8.3s |
| remove-every-chart-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 8.4s |
| change-color-column | atomic-do | pass | pass | pass | pass | **PASS** | 8.6s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- `wrong-tool-trap-fromtype-direct`: passed rubric but emitted no exec blocks (expected_path may not have required one).
