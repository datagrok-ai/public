# Eval report — selection

Started: 2026-05-18T12:28:51.099Z
Duration: 34.4s
Model: `sonnet`
System prompt size: 27,124 chars
Inlined skills: `datagrok-exec`, `selection`

## Summary

- **Total prompts:** 18
- **Passed:** 18
- **Failed:** 0
- **Pass rate:** 100.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 5 | 5 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 4 | 4 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 34.4s wall
- Mean: 6.6s · Median: 5.9s · p95: 12.3s
- Slowest 3: `describe-current-selection` (12.3s), `reset-select-all-molecule-column` (11.5s), `select-age-gt-40` (8.9s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`selection`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| select-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 8.9s |
| select-rows-0-10 | atomic-do | pass | pass | pass | pass | **PASS** | 6.8s |
| select-top-5-by-activity | atomic-do | pass | pass | pass | pass | **PASS** | 7.1s |
| select-ic50-under-100 | atomic-do | pass | pass | pass | pass | **PASS** | 8.8s |
| select-all-rows | atomic-do | pass | pass | pass | pass | **PASS** | 4.3s |
| clear-the-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 5.1s |
| deselect-everything | atomic-undo | pass | pass | pass | pass | **PASS** | 5.0s |
| invert-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 4.6s |
| remove-rows-5-and-7 | atomic-undo | pass | pass | pass | pass | **PASS** | 5.1s |
| show-selected-as-new-table | combine | pass | pass | pass | pass | **PASS** | 5.3s |
| select-all-visible-rows | combine | pass | pass | pass | pass | **PASS** | 5.4s |
| filter-to-selected-rows | combine | pass | pass | pass | pass | **PASS** | 7.3s |
| wrong-tool-trap-clear-via-setall-true | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.8s |
| wrong-tool-trap-select-all-not-clear | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.6s |
| wrong-tool-trap-highlight-row-5-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.9s |
| reset-deselect-and-default | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 6.8s |
| reset-select-all-molecule-column | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 11.5s |
| describe-current-selection | combine | pass | pass | pass | pass | **PASS** | 12.3s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
