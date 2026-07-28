# Eval report — df-and-columns

Started: 2026-05-18T11:42:39.968Z
Duration: 36.4s
Model: `sonnet`
System prompt size: 36,054 chars
Inlined skills: `datagrok-exec`, `df-and-columns`

## Summary

- **Total prompts:** 18
- **Passed:** 18
- **Failed:** 0
- **Pass rate:** 100.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 6 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 36.4s wall
- Mean: 7.2s · Median: 7.1s · p95: 17.7s
- Slowest 3: `describe-activity-col` (17.7s), `describe-and-add` (13.1s), `add-logmw-from-mw` (10.8s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`df-and-columns`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| find-molecule-col | atomic-do | pass | pass | pass | pass | **PASS** | 8.8s |
| describe-activity-col | atomic-do | pass | pass | pass | pass | **PASS** | 17.7s |
| add-logmw-from-mw | atomic-do | pass | pass | pass | pass | **PASS** | 10.8s |
| add-empty-float-with-units | atomic-do | pass | pass | pass | pass | **PASS** | 7.4s |
| set-friendly-name | atomic-do | pass | pass | pass | pass | **PASS** | 4.5s |
| color-coding-linear | atomic-do | pass | pass | pass | pass | **PASS** | 4.8s |
| remove-just-added | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| clear-friendly-name | atomic-undo | pass | pass | pass | pass | **PASS** | 4.5s |
| rename-back | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| color-coding-off | atomic-undo | pass | pass | pass | pass | **PASS** | 8.6s |
| describe-and-add | combine | pass | pass | pass | pass | **PASS** | 13.1s |
| clone-filtered-subset | combine | pass | pass | pass | pass | **PASS** | 4.8s |
| find-and-describe-numeric | combine | pass | pass | pass | pass | **PASS** | 6.6s |
| wrong-tool-trap-semtype | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 7.1s |
| wrong-tool-trap-float-literal | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.5s |
| wrong-tool-trap-color-tag | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.8s |
| reset-all-meta | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 7.4s |
| remove-session-calculated | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 7.4s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
