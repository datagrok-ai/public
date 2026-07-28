# Eval report — viewers

Started: 2026-05-18T11:35:14.693Z
Duration: 36.4s
Model: `sonnet`
System prompt size: 33,235 chars
Inlined skills: `datagrok-exec`, `viewers`

## Summary

- **Total prompts:** 20
- **Passed:** 17
- **Failed:** 3
- **Pass rate:** 85.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 8 | 75% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 2 | 3 | 67% |

### Latency

- Total: 36.4s wall
- Mean: 6.9s · Median: 5.3s · p95: 18.7s
- Slowest 3: `add-bar-chart-count-by-category` (18.7s), `wrong-tool-trap-chem-space` (15.2s), `add-histogram-split-by-class` (11.3s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`viewers`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| add-scatter-height-weight | atomic-do | pass | pass | pass | pass | **PASS** | 4.3s |
| add-histogram-age | atomic-do | pass | pass | pass | pass | **PASS** | 4.7s |
| add-line-chart-date-total | atomic-do | pass | pass | pass | pass | **PASS** | 6.6s |
| add-histogram-split-by-class | atomic-do | pass | fail | fail | fail | **FAIL** | 11.3s |
| add-bar-chart-count-by-category | atomic-do | pass | pass | fail | pass | **FAIL** | 18.7s |
| add-box-plot-activity-by-class | atomic-do | pass | pass | pass | pass | **PASS** | 7.2s |
| add-correlation-plot-numeric | atomic-do | pass | pass | pass | pass | **PASS** | 6.1s |
| close-scatter-plot | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| close-all-except-grid | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| remove-histogram | atomic-undo | pass | pass | pass | pass | **PASS** | 5.1s |
| close-every-chart | atomic-undo | pass | pass | pass | pass | **PASS** | 4.4s |
| scatter-mw-logp-colored-regression | combine | pass | pass | pass | pass | **PASS** | 5.2s |
| add-then-close-sequenced | combine | pass | pass | pass | pass | **PASS** | 4.8s |
| replace-scatter-with-heatmap | combine | pass | pass | pass | pass | **PASS** | 8.2s |
| wrong-tool-trap-chem-space | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 15.2s |
| wrong-tool-trap-deprecated-shorthand | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.5s |
| wrong-tool-trap-fromtype-direct | wrong-tool-trap | fail | fail | fail | pass | **FAIL** | 8.6s |
| reset-to-just-grid | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.4s |
| remove-every-chart-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 5.3s |
| change-color-column | atomic-do | pass | pass | pass | pass | **PASS** | 4.7s |

## Failed prompts (deep)

### add-histogram-split-by-class — atomic-do

**Prompt:** Show me the distribution of activities per compound class
**Intent:** Histogram of 'activity' split by category column. splitColumnName: 'compound_class', valueColumnName: 'activity'. Box plot would also visualize this but the user said 'distribution' which is histogram-canonical.

**Rubric:**

- path: PASS
- helpers: FAIL (`view\.addViewer\(\s*DG\.VIEWER\.HISTOGRAM`)
- symbols: FAIL (`DG\.VIEWER\.HISTOGRAM`, `splitColumnName`)
- forbidden: FAIL (`DG\.VIEWER\.BOX_PLOT`)

**Hypothesis:** Claude bypassed the expected helper(s) (`view\.addViewer\(\s*DG\.VIEWER\.HISTOGRAM`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`DG\.VIEWER\.HISTOGRAM`, `splitColumnName`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill. Forbidden pattern(s) matched (`DG\.VIEWER\.BOX_PLOT`). Skill's anti-pattern callout isn't strong enough — consider making the wrong path throw, or surface a more obvious "DO NOT" callout near the relevant example.

**Emitted exec blocks:**

```javascript
// Return column names and types to identify activity and class columns
const cols = Array.from(t.columns).map(c => ({name: c.name, type: c.type}));
return ui.tableFromMap(Object.fromEntries(cols.map(c => [c.name, c.type])));
```

```javascript
view.addViewer(DG.VIEWER.BOX_PLOT, {
  valueColumnName: 'activity',      // ← replace with your activity column
  categoryColumnName: 'class',      // ← replace with your compound class column
});
```

<details><summary>Full response text</summary>

```
I'll inspect the current table to find the right columns automatically.

```datagrok-exec
// Return column names and types to identify activity and class columns
const cols = Array.from(t.columns).map(c => ({name: c.name, type: c.type}));
return ui.tableFromMap(Object.fromEntries(cols.map(c => [c.name, c.type])));
```

Could you let me know which columns represent **activity** and **compound class** in your table? Once I have those, I'll add a Box Plot (best for showing distributions per category) to visualize the activity distribution per compound class.

Alternatively, if you can share the column names, I can add the viewer right away — for example:

```datagrok-exec
view.addViewer(DG.VIEWER.BOX_PLOT, {
  valueColumnName: 'activity',      // ← replace with your activity column
  categoryColumnName: 'class',      // ← replace with your compound class column
});
```

Please confirm the column names and I'll render it immediately!
```

</details>

### add-bar-chart-count-by-category — atomic-do

**Prompt:** Add a bar chart of count by category
**Intent:** Bar chart aggregating count by a category column. splitColumnName: 'category' is the canonical shape; valueAggrType: 'count' commonly accompanies it.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`splitColumnName\s*:\s*['"]category['"]`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`splitColumnName\s*:\s*['"]category['"]`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
view.addViewer(DG.VIEWER.BAR_CHART, {
  valueAggrType: 'count',
});
```

<details><summary>Full response text</summary>

```
I'll add a bar chart using the default categorical column for grouping. Datagrok will automatically pick the first suitable column.

```datagrok-exec
view.addViewer(DG.VIEWER.BAR_CHART, {
  valueAggrType: 'count',
});
```

A **Bar Chart** has been added to your view, showing **count by category**. Datagrok will automatically use the first categorical column as the split (category) axis.

If you'd like to change which column is used for grouping, just let me know the column name and I'll update it!
```

</details>

### wrong-tool-trap-fromtype-direct — wrong-tool-trap

**Prompt:** Use DG.Viewer.fromType for new viewers
**Intent:** Trap: DG.Viewer.fromType creates a DETACHED viewer that the caller must dock. On a TableView, view.addViewer is correct — it creates + docks in one call.
**Trap:** DG.Viewer.fromType is only correct for off-view rendering (e.g. embedding in a returned HTMLElement). On a TableView, use view.addViewer(type, opts).

**Rubric:**

- path: FAIL (`>=1 datagrok-exec block`)
- helpers: FAIL (`view\.addViewer`)
- symbols: FAIL (`view\.addViewer`)
- forbidden: PASS

**Hypothesis:** Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class. Claude bypassed the expected helper(s) (`view\.addViewer`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`view\.addViewer`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:** (none)

<details><summary>Full response text</summary>

```
Noted! I'll use `DG.Viewer.fromType(type, df, opts?)` instead of `view.addViewer(...)` when creating new viewers in `datagrok-exec` blocks.

Just to confirm the distinction:

- **`DG.Viewer.fromType(type, df, opts?)`** — creates a standalone (detached) viewer; you get back the viewer object and are responsible for placing/returning its `.root`.
- **`view.addViewer(type, opts?)`** — creates **and docks** the viewer on the current `TableView` in one call; no placement needed.

Which should I use by default depends on the goal:

| Goal | Preferred call |
|---|---|
| Add a viewer to the current view (most common) | `view.addViewer(type, opts)` |
| Embed a viewer in a custom UI / return it to chat | `DG.Viewer.fromType(type, df, opts)` → return `viewer.root` |

Do you want me to **always** use `DG.Viewer.fromType` (even when docking to the view), or only when creating standalone/embedded viewers?
```

</details>

## Suggested skill improvements

- Expected helper `view\.addViewer\(\s*DG\.VIEWER\.HISTOGRAM` not reached for: add-histogram-split-by-class. Consider tightening the skill's routing table or example density for this case.
- Expected helper `view\.addViewer` not reached for: wrong-tool-trap-fromtype-direct. Consider tightening the skill's routing table or example density for this case.
- Forbidden pattern `DG\.VIEWER\.BOX_PLOT` slipped through for: add-histogram-split-by-class. Anti-pattern callout needs to be more prominent — consider a "WRONG / RIGHT" pair right before the relevant example.

## Anomalies

- None.
