# Eval report — selection

Started: 2026-05-14T10:43:38.233Z
Duration: 33.9s
Model: `sonnet`
System prompt size: 33,974 chars
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

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`selection`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| select-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 5.4s |
| select-rows-0-10 | atomic-do | pass | pass | pass | pass | **PASS** | 4.2s |
| select-top-5-by-activity | atomic-do | pass | pass | pass | pass | **PASS** | 12.2s |
| select-ic50-under-100 | atomic-do | pass | pass | pass | pass | **PASS** | 10.9s |
| select-all-rows | atomic-do | pass | pass | pass | pass | **PASS** | 7.9s |
| clear-the-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 3.7s |
| deselect-everything | atomic-undo | pass | pass | pass | pass | **PASS** | 3.7s |
| invert-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| remove-rows-5-and-7 | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| show-selected-as-new-table | combine | pass | pass | pass | pass | **PASS** | 4.3s |
| select-all-visible-rows | combine | pass | pass | pass | pass | **PASS** | 3.7s |
| filter-to-selected-rows | combine | pass | pass | pass | pass | **PASS** | 4.1s |
| wrong-tool-trap-clear-via-setall-true | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.8s |
| wrong-tool-trap-select-all-not-clear | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.3s |
| wrong-tool-trap-highlight-row-5-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.4s |
| reset-deselect-and-default | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 14.8s |
| reset-select-all-molecule-column | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 10.5s |
| describe-current-selection | combine | pass | pass | pass | pass | **PASS** | 3.5s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
