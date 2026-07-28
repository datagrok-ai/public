# Eval report — df-and-columns

Started: 2026-05-18T12:28:07.063Z
Duration: 54.1s
Model: `sonnet`
System prompt size: 30,074 chars
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

- Total: 54.1s wall
- Mean: 9.6s · Median: 6.5s · p95: 23.9s
- Slowest 3: `describe-activity-col` (23.9s), `find-and-describe-numeric` (21.5s), `describe-and-add` (19.7s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`df-and-columns`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| find-molecule-col | atomic-do | pass | pass | pass | pass | **PASS** | 9.3s |
| describe-activity-col | atomic-do | pass | pass | pass | pass | **PASS** | 23.9s |
| add-logmw-from-mw | atomic-do | pass | pass | pass | pass | **PASS** | 9.1s |
| add-empty-float-with-units | atomic-do | pass | pass | pass | pass | **PASS** | 4.6s |
| set-friendly-name | atomic-do | pass | pass | pass | pass | **PASS** | 4.1s |
| color-coding-linear | atomic-do | pass | pass | pass | pass | **PASS** | 7.0s |
| remove-just-added | atomic-undo | pass | pass | pass | pass | **PASS** | 4.7s |
| clear-friendly-name | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| rename-back | atomic-undo | pass | pass | pass | pass | **PASS** | 4.9s |
| color-coding-off | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| describe-and-add | combine | pass | pass | pass | pass | **PASS** | 19.7s |
| clone-filtered-subset | combine | pass | pass | pass | pass | **PASS** | 18.3s |
| find-and-describe-numeric | combine | pass | pass | pass | pass | **PASS** | 21.5s |
| wrong-tool-trap-semtype | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.1s |
| wrong-tool-trap-float-literal | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.5s |
| wrong-tool-trap-color-tag | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.7s |
| reset-all-meta | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 3.8s |
| remove-session-calculated | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 17.5s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
