# Eval report — grid-customization

Started: 2026-05-14T11:18:24.146Z
Duration: 24.6s
Model: `sonnet`
System prompt size: 34,180 chars
Inlined skills: `datagrok-exec`, `grid-customization`

## Summary

- **Total prompts:** 19
- **Passed:** 18
- **Failed:** 1
- **Pass rate:** 94.7%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 7 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 2 | 3 | 67% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`grid-customization`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| sort-activity-desc | atomic-do | pass | pass | pass | pass | **PASS** | 4.6s |
| hide-index-column | atomic-do | pass | pass | pass | pass | **PASS** | 4.6s |
| pin-smiles-left | atomic-do | pass | pass | pass | pass | **PASS** | 6.2s |
| widen-molecule-column | atomic-do | pass | pass | pass | pass | **PASS** | 5.2s |
| color-activity-red-green | atomic-do | pass | pass | pass | pass | **PASS** | 6.1s |
| format-ic50-two-decimals | atomic-do | pass | pass | pass | pass | **PASS** | 4.0s |
| show-only-five-columns | atomic-do | pass | pass | pass | pass | **PASS** | 4.9s |
| clear-sort | atomic-undo | pass | pass | pass | pass | **PASS** | 3.6s |
| show-all-hidden | atomic-undo | pass | pass | pass | pass | **PASS** | 5.2s |
| unpin-smiles | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| color-off-activity | atomic-undo | pass | pass | pass | pass | **PASS** | 4.5s |
| hit-list-multi-step | combine | pass | pass | pass | pass | **PASS** | 7.2s |
| reset-grid-keep-viewers | combine | pass | pass | pass | pass | **PASS** | 4.7s |
| sort-class-then-activity | combine | pass | pass | pass | pass | **PASS** | 4.4s |
| trap-sortbycolumns | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.7s |
| trap-gridcol-name-rename | wrong-tool-trap | pass | fail | pass | pass | **FAIL** | 6.7s |
| trap-reset-layout-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.2s |
| reset-all-customization | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 5.0s |
| clear-everything-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.2s |

## Failed prompts (deep)

### trap-gridcol-name-rename — wrong-tool-trap

**Prompt:** Rename the activity column in the grid display to 'Potency'
**Intent:** User asks for grid-display rename. gridCol.name = 'Potency' actually renames the underlying DF column (destructive — breaks every reference). Correct: col.meta.friendlyName = 'Potency' or setColumnMeta(col, {friendlyName: 'Potency'}). The skill explicitly traps this.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.setColumnMeta`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.setColumnMeta`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const col = t.col('activity');
if (!col) throw new Error("Column 'activity' not found in the current table.");
col.meta.friendlyName = 'Potency';
```

<details><summary>Full response text</summary>

```

```datagrok-exec
const col = t.col('activity');
if (!col) throw new Error("Column 'activity' not found in the current table.");
col.meta.friendlyName = 'Potency';
```

Done! The **activity** column is now displayed as **Potency** in the grid (and all other viewers). The underlying column name remains `activity`, so any existing formulas or scripts referencing it will continue to work.
```

</details>

## Suggested skill improvements

- Expected helper `grokky.setColumnMeta` not reached for: trap-gridcol-name-rename. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
