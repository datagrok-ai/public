# Eval report — filtering

Started: 2026-05-18T13:12:48.858Z
Duration: 74.7s
Model: `sonnet`
System prompt size: 27,949 chars
Inlined skills: `datagrok-exec`, `filtering`

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

- Total: 74.7s wall
- Mean: 11.2s · Median: 6.4s · p95: 65.0s
- Slowest 3: `filter-smarts-pattern` (65.0s), `filter-mw-lt-500` (19.6s), `cross-skill-selected-then-describe` (16.3s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`filtering`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| filter-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 11.3s |
| filter-mw-lt-500 | atomic-do | pass | pass | pass | pass | **PASS** | 19.6s |
| filter-substructure-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 6.4s |
| filter-categorical-in-set | atomic-do | pass | pass | pass | pass | **PASS** | 9.7s |
| filter-name-contains | atomic-do | pass | pass | pass | pass | **PASS** | 4.5s |
| filter-smarts-pattern | atomic-do | pass | pass | pass | pass | **PASS** | 65.0s |
| clear-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 4.7s |
| show-all-rows | atomic-undo | pass | pass | pass | pass | **PASS** | 4.2s |
| remove-substructure-only | atomic-undo | pass | pass | pass | pass | **PASS** | 8.2s |
| invert-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 4.2s |
| filter-mw-and-logp | combine | pass | pass | pass | pass | **PASS** | 5.1s |
| filtered-as-new-table | combine | pass | pass | pass | pass | **PASS** | 4.5s |
| cross-skill-selected-then-describe | combine | pass | pass | pass | pass | **PASS** | 16.3s |
| wrong-tool-trap-remove-as-delete | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.6s |
| wrong-tool-trap-smiles-as-molblock | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 7.3s |
| wrong-tool-trap-hide-vs-drop | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.6s |
| reset-clear-all-filters-and-selections | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.0s |
| go-back-to-seeing-every-row | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 16.2s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
