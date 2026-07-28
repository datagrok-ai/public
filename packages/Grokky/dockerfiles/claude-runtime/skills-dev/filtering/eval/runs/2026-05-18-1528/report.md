# Eval report — filtering

Started: 2026-05-18T12:28:10.645Z
Duration: 54.1s
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

- Total: 54.1s wall
- Mean: 8.0s · Median: 4.8s · p95: 37.0s
- Slowest 3: `cross-skill-selected-then-describe` (37.0s), `wrong-tool-trap-smiles-as-molblock` (16.4s), `filter-smarts-pattern` (11.5s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`filtering`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| filter-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 5.6s |
| filter-mw-lt-500 | atomic-do | pass | pass | pass | pass | **PASS** | 8.8s |
| filter-substructure-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 7.2s |
| filter-categorical-in-set | atomic-do | pass | pass | pass | pass | **PASS** | 4.8s |
| filter-name-contains | atomic-do | pass | pass | pass | pass | **PASS** | 4.2s |
| filter-smarts-pattern | atomic-do | pass | pass | pass | pass | **PASS** | 11.5s |
| clear-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 4.6s |
| show-all-rows | atomic-undo | pass | pass | pass | pass | **PASS** | 4.6s |
| remove-substructure-only | atomic-undo | pass | pass | pass | pass | **PASS** | 8.0s |
| invert-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 3.7s |
| filter-mw-and-logp | combine | pass | pass | pass | pass | **PASS** | 4.7s |
| filtered-as-new-table | combine | pass | pass | pass | pass | **PASS** | 4.7s |
| cross-skill-selected-then-describe | combine | pass | pass | pass | pass | **PASS** | 37.0s |
| wrong-tool-trap-remove-as-delete | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.1s |
| wrong-tool-trap-smiles-as-molblock | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 16.4s |
| wrong-tool-trap-hide-vs-drop | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.2s |
| reset-clear-all-filters-and-selections | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.5s |
| go-back-to-seeing-every-row | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 3.7s |

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
