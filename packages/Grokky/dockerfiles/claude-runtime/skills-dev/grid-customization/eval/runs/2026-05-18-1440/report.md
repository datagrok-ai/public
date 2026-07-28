# Eval report — grid-customization

Started: 2026-05-18T11:40:31.893Z
Duration: 31.1s
Model: `sonnet`
System prompt size: 25,838 chars
Inlined skills: `datagrok-exec`, `grid-customization`

## Summary

- **Total prompts:** 19
- **Passed:** 18
- **Failed:** 1
- **Pass rate:** 94.7%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 7 | 86% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 31.1s wall
- Mean: 5.7s · Median: 4.6s · p95: 15.0s
- Slowest 3: `sort-activity-desc` (15.0s), `widen-molecule-column` (10.7s), `trap-sortbycolumns` (8.1s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`grid-customization`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| sort-activity-desc | atomic-do | pass | pass | pass | pass | **PASS** | 15.0s |
| hide-index-column | atomic-do | pass | pass | fail | pass | **FAIL** | 5.3s |
| pin-smiles-left | atomic-do | pass | pass | pass | pass | **PASS** | 3.9s |
| widen-molecule-column | atomic-do | pass | pass | pass | pass | **PASS** | 10.7s |
| color-activity-red-green | atomic-do | pass | pass | pass | pass | **PASS** | 4.8s |
| format-ic50-two-decimals | atomic-do | pass | pass | pass | pass | **PASS** | 4.3s |
| show-only-five-columns | atomic-do | pass | pass | pass | pass | **PASS** | 4.0s |
| clear-sort | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| show-all-hidden | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| unpin-smiles | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| color-off-activity | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| hit-list-multi-step | combine | pass | pass | pass | pass | **PASS** | 4.6s |
| reset-grid-keep-viewers | combine | pass | pass | pass | pass | **PASS** | 5.5s |
| sort-class-then-activity | combine | pass | pass | pass | pass | **PASS** | 5.0s |
| trap-sortbycolumns | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 8.1s |
| trap-gridcol-name-rename | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.1s |
| trap-reset-layout-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.0s |
| reset-all-customization | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 5.7s |
| clear-everything-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 7.6s |

## Failed prompts (deep)

### hide-index-column — atomic-do

**Prompt:** Hide the index column
**Intent:** Display-only hide via gridCol.visible = false or grid.columns.setVisible([...]) (excluding 'index'). Must not remove the column from the DataFrame. Datagrok's row-header column is named '#' — Claude may correctly target that instead of 'index'.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`index|['"]#['"]`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`index|['"]#['"]`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
view.grid.col('').visible = false;
```

<details><summary>Full response text</summary>

```

```datagrok-exec
view.grid.col('').visible = false;
```

The row-index column (the unlabeled first column) is now hidden.
```

</details>

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
