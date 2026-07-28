# Eval report — grid-customization

Started: 2026-05-14T11:15:55.066Z
Duration: 27.5s
Model: `sonnet`
System prompt size: 34,180 chars
Inlined skills: `datagrok-exec`, `grid-customization`

## Summary

- **Total prompts:** 19
- **Passed:** 11
- **Failed:** 8
- **Pass rate:** 57.9%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 4 | 7 | 57% |
| atomic-undo | 2 | 4 | 50% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 1 | 2 | 50% |
| wrong-tool-trap | 1 | 3 | 33% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`grid-customization`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| sort-activity-desc | atomic-do | pass | pass | pass | pass | **PASS** | 4.9s |
| hide-index-column | atomic-do | pass | pass | fail | pass | **FAIL** | 4.6s |
| pin-smiles-left | atomic-do | pass | fail | pass | pass | **FAIL** | 4.3s |
| widen-molecule-column | atomic-do | pass | pass | pass | pass | **PASS** | 4.5s |
| color-activity-red-green | atomic-do | pass | pass | pass | pass | **PASS** | 5.9s |
| format-ic50-two-decimals | atomic-do | pass | fail | pass | pass | **FAIL** | 4.8s |
| show-only-five-columns | atomic-do | pass | pass | pass | pass | **PASS** | 4.9s |
| clear-sort | atomic-undo | pass | pass | pass | pass | **PASS** | 4.5s |
| show-all-hidden | atomic-undo | pass | fail | pass | pass | **FAIL** | 6.7s |
| unpin-smiles | atomic-undo | pass | fail | pass | pass | **FAIL** | 4.7s |
| color-off-activity | atomic-undo | pass | pass | pass | pass | **PASS** | 4.6s |
| hit-list-multi-step | combine | pass | pass | pass | pass | **PASS** | 5.9s |
| reset-grid-keep-viewers | combine | pass | pass | pass | pass | **PASS** | 4.5s |
| sort-class-then-activity | combine | pass | pass | pass | pass | **PASS** | 4.3s |
| trap-sortbycolumns | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.5s |
| trap-gridcol-name-rename | wrong-tool-trap | pass | fail | pass | pass | **FAIL** | 5.5s |
| trap-reset-layout-ambiguous | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 3.7s |
| reset-all-customization | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.4s |
| clear-everything-fresh-start | reset-or-cleanup | pass | fail | pass | pass | **FAIL** | 6.1s |

## Failed prompts (deep)

### hide-index-column — atomic-do

**Prompt:** Hide the index column
**Intent:** Display-only hide via configureGrid({hide: ['index']}) or gridCol.visible = false. Must not remove the column from the DataFrame.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`index`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`index`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
grokky.configureGrid(view, {hide: ['#']});
```

<details><summary>Full response text</summary>

```
```datagrok-exec
grokky.configureGrid(view, {hide: ['#']});
```
```

</details>

### pin-smiles-left — atomic-do

**Prompt:** Pin the SMILES column to the left
**Intent:** Built-in pin via gridCol.pin() (left-only) — either grokky.pinColumn(view, 'SMILES') or configureGrid({pin: ['SMILES']}). PowerGrid plugin path should not be required.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.configureGrid`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.configureGrid`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
grokky.pinColumn(view, 'SMILES');
```

<details><summary>Full response text</summary>

```
```datagrok-exec
grokky.pinColumn(view, 'SMILES');
```

The **SMILES** column has been pinned to the left of the grid.
```

</details>

### format-ic50-two-decimals — atomic-do

**Prompt:** Set the format of IC50 to 2 decimal places
**Intent:** configureGrid({formats: {IC50: '0.00'}}) is canonical for this skill. setColumnMeta from df-and-columns also works and is preferred when in a metadata pipeline. Either passes.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.setColumnMeta`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.setColumnMeta`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
grokky.configureGrid(view, {formats: {IC50: '0.00'}});
```

<details><summary>Full response text</summary>

```
```datagrok-exec
grokky.configureGrid(view, {formats: {IC50: '0.00'}});
```

Done! The **IC50** column is now formatted to 2 decimal places.
```

</details>

### show-all-hidden — atomic-undo

**Prompt:** Show all hidden columns
**Intent:** Restore visibility for every column. grokky.resetGrid(view, {visibility: true, widths: false, sort: false, colors: false}) is the precise tool; a loop over grid.columns toggling visible=true is also acceptable. resetLayout is too blunt.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.configureGrid`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.configureGrid`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
grokky.resetGrid(view, {visibility: true, widths: false, sort: false, colors: false});
```

<details><summary>Full response text</summary>

```

```datagrok-exec
grokky.resetGrid(view, {visibility: true, widths: false, sort: false, colors: false});
```

All hidden columns are now visible again.
```

</details>

### unpin-smiles — atomic-undo

**Prompt:** Unpin the SMILES column
**Intent:** grokky.unpinColumn(view, 'SMILES') or configureGrid({unpin: ['SMILES']}) or gridCol.unpin() directly.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.configureGrid`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.configureGrid`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
grokky.unpinColumn(view, 'SMILES');
```

<details><summary>Full response text</summary>

```
```datagrok-exec
grokky.unpinColumn(view, 'SMILES');
```

The **SMILES** column has been unpinned.
```

</details>

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

### trap-reset-layout-ambiguous — wrong-tool-trap

**Prompt:** Reset the layout to just the grid
**Intent:** Ambiguous: 'reset the layout to just the grid' literally maps to view.resetLayout() (which closes viewers AND resets the grid), but a user might mean grokky.resetGrid alone. Pass if Claude either (a) picks view.resetLayout() as the natural interpretation, (b) calls grokky.resetGrid + closeAllViewers, or (c) surfaces the distinction in prose before picking.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.resetGrid`, `grokky.closeAllViewers`)
- symbols: FAIL (`grid|viewer`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.resetGrid`, `grokky.closeAllViewers`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`grid|viewer`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
view.resetLayout();
```

<details><summary>Full response text</summary>

```
```datagrok-exec
view.resetLayout();
```

The layout has been reset — all additional viewers have been closed and the view now shows just the grid.
```

</details>

### clear-everything-fresh-start — reset-or-cleanup

**Prompt:** Clear everything and start fresh
**Intent:** Broader than reset-all-customization — 'everything' arguably includes added viewers. Pass if grokky.resetGrid is called; closing other viewers is a bonus. view.resetLayout() also passes (it's the platform's blunt instrument for exactly this). Only fail if Claude does nothing or only changes one aspect.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.closeAllViewers`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.closeAllViewers`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
grokky.resetGrid(view);
```

<details><summary>Full response text</summary>

```
```datagrok-exec
grokky.resetGrid(view);
```

Done! The grid has been reset — visibility, widths, sort, and color coding are all cleared back to defaults.
```

</details>

## Suggested skill improvements

- Expected helper `grokky.configureGrid` not reached for: pin-smiles-left, show-all-hidden, unpin-smiles. Consider tightening the skill's routing table or example density for this case.
- Expected helper `grokky.setColumnMeta` not reached for: format-ic50-two-decimals, trap-gridcol-name-rename. Consider tightening the skill's routing table or example density for this case.
- Expected helper `grokky.resetGrid` not reached for: trap-reset-layout-ambiguous. Consider tightening the skill's routing table or example density for this case.
- Expected helper `grokky.closeAllViewers` not reached for: trap-reset-layout-ambiguous, clear-everything-fresh-start. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
