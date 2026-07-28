# Eval report — df-and-columns

Started: 2026-05-14T10:03:51.454Z
Duration: 38.5s
Model: `sonnet`
System prompt size: 32,677 chars
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

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`df-and-columns`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| find-molecule-col | atomic-do | pass | pass | pass | pass | **PASS** | 13.7s |
| describe-activity-col | atomic-do | pass | pass | pass | pass | **PASS** | 16.9s |
| add-logmw-from-mw | atomic-do | pass | pass | pass | pass | **PASS** | 7.5s |
| add-empty-float-with-units | atomic-do | pass | pass | pass | pass | **PASS** | 5.0s |
| set-friendly-name | atomic-do | pass | pass | pass | pass | **PASS** | 4.2s |
| color-coding-linear | atomic-do | pass | pass | pass | pass | **PASS** | 5.0s |
| remove-just-added | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| clear-friendly-name | atomic-undo | pass | pass | pass | pass | **PASS** | 7.4s |
| rename-back | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| color-coding-off | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| describe-and-add | combine | pass | pass | pass | pass | **PASS** | 13.1s |
| clone-filtered-subset | combine | pass | pass | pass | pass | **PASS** | 11.5s |
| find-and-describe-numeric | combine | pass | pass | pass | pass | **PASS** | 4.7s |
| wrong-tool-trap-semtype | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 11.7s |
| wrong-tool-trap-float-literal | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.8s |
| wrong-tool-trap-color-tag | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.9s |
| reset-all-meta | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 9.6s |
| remove-session-calculated | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.7s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
