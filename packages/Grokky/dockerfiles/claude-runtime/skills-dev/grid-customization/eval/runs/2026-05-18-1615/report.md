# Eval report — grid-customization

Started: 2026-05-18T13:15:23.182Z
Duration: 48.0s
Model: `sonnet`
System prompt size: 22,721 chars
Inlined skills: `datagrok-exec`, `grid-customization`

## Summary

- **Total prompts:** 19
- **Passed:** 19
- **Failed:** 0
- **Pass rate:** 100.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 7 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 48.0s wall
- Mean: 9.1s · Median: 9.2s · p95: 18.3s
- Slowest 3: `widen-molecule-column` (18.3s), `hit-list-multi-step` (16.3s), `show-all-hidden` (15.9s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`grid-customization`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| sort-activity-desc | atomic-do | pass | pass | pass | pass | **PASS** | 4.3s |
| hide-index-column | atomic-do | pass | pass | pass | pass | **PASS** | 15.2s |
| pin-smiles-left | atomic-do | pass | pass | pass | pass | **PASS** | 14.5s |
| widen-molecule-column | atomic-do | pass | pass | pass | pass | **PASS** | 18.3s |
| color-activity-red-green | atomic-do | pass | pass | pass | pass | **PASS** | 9.5s |
| format-ic50-two-decimals | atomic-do | pass | pass | pass | pass | **PASS** | 7.2s |
| show-only-five-columns | atomic-do | pass | pass | pass | pass | **PASS** | 3.9s |
| clear-sort | atomic-undo | pass | pass | pass | pass | **PASS** | 9.5s |
| show-all-hidden | atomic-undo | pass | pass | pass | pass | **PASS** | 15.9s |
| unpin-smiles | atomic-undo | pass | pass | pass | pass | **PASS** | 3.8s |
| color-off-activity | atomic-undo | pass | pass | pass | pass | **PASS** | 3.6s |
| hit-list-multi-step | combine | pass | pass | pass | pass | **PASS** | 16.3s |
| reset-grid-keep-viewers | combine | pass | pass | pass | pass | **PASS** | 6.0s |
| sort-class-then-activity | combine | pass | pass | pass | pass | **PASS** | 4.1s |
| trap-sortbycolumns | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 9.7s |
| trap-gridcol-name-rename | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 9.2s |
| trap-reset-layout-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.9s |
| reset-all-customization | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 9.5s |
| clear-everything-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 6.9s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
