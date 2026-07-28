# Eval report — viewers

Started: 2026-05-18T12:29:15.264Z
Duration: 45.4s
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

- Total: 45.4s wall
- Mean: 8.4s · Median: 7.7s · p95: 18.6s
- Slowest 3: `replace-scatter-with-heatmap` (18.6s), `add-histogram-split-by-class` (13.9s), `wrong-tool-trap-chem-space` (13.1s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`viewers`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| add-scatter-height-weight | atomic-do | pass | pass | pass | pass | **PASS** | 7.3s |
| add-histogram-age | atomic-do | pass | pass | pass | pass | **PASS** | 5.5s |
| add-line-chart-date-total | atomic-do | pass | pass | pass | pass | **PASS** | 6.5s |
| add-histogram-split-by-class | atomic-do | pass | pass | pass | pass | **PASS** | 13.9s |
| add-bar-chart-count-by-category | atomic-do | pass | pass | pass | pass | **PASS** | 5.1s |
| add-box-plot-activity-by-class | atomic-do | pass | pass | pass | pass | **PASS** | 5.9s |
| add-correlation-plot-numeric | atomic-do | pass | pass | pass | pass | **PASS** | 6.8s |
| close-scatter-plot | atomic-undo | pass | pass | pass | pass | **PASS** | 5.7s |
| close-all-except-grid | atomic-undo | pass | pass | pass | pass | **PASS** | 9.9s |
| remove-histogram | atomic-undo | pass | pass | pass | pass | **PASS** | 10.5s |
| close-every-chart | atomic-undo | pass | pass | pass | pass | **PASS** | 8.1s |
| scatter-mw-logp-colored-regression | combine | pass | pass | pass | pass | **PASS** | 7.7s |
| add-then-close-sequenced | combine | pass | pass | pass | pass | **PASS** | 5.7s |
| replace-scatter-with-heatmap | combine | pass | pass | pass | pass | **PASS** | 18.6s |
| wrong-tool-trap-chem-space | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 13.1s |
| wrong-tool-trap-deprecated-shorthand | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.7s |
| wrong-tool-trap-fromtype-direct | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 9.1s |
| reset-to-just-grid | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 3.8s |
| remove-every-chart-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 11.7s |
| change-color-column | atomic-do | pass | pass | pass | pass | **PASS** | 8.5s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- `wrong-tool-trap-fromtype-direct`: passed rubric but emitted no exec blocks (expected_path may not have required one).
