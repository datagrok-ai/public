# Eval report — viewers

Started: 2026-05-18T11:40:27.755Z
Duration: 39.0s
Model: `sonnet`
System prompt size: 33,235 chars
Inlined skills: `datagrok-exec`, `viewers`

## Summary

- **Total prompts:** 20
- **Passed:** 19
- **Failed:** 1
- **Pass rate:** 95.0%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 8 | 88% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 39.0s wall
- Mean: 6.7s · Median: 5.5s · p95: 14.9s
- Slowest 3: `wrong-tool-trap-chem-space` (14.9s), `close-all-except-grid` (10.9s), `add-box-plot-activity-by-class` (10.5s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`viewers`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| add-scatter-height-weight | atomic-do | pass | pass | pass | pass | **PASS** | 4.2s |
| add-histogram-age | atomic-do | pass | pass | pass | pass | **PASS** | 5.3s |
| add-line-chart-date-total | atomic-do | pass | pass | pass | pass | **PASS** | 4.5s |
| add-histogram-split-by-class | atomic-do | pass | fail | fail | pass | **FAIL** | 10.1s |
| add-bar-chart-count-by-category | atomic-do | pass | pass | pass | pass | **PASS** | 5.2s |
| add-box-plot-activity-by-class | atomic-do | pass | pass | pass | pass | **PASS** | 10.5s |
| add-correlation-plot-numeric | atomic-do | pass | pass | pass | pass | **PASS** | 5.6s |
| close-scatter-plot | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| close-all-except-grid | atomic-undo | pass | pass | pass | pass | **PASS** | 10.9s |
| remove-histogram | atomic-undo | pass | pass | pass | pass | **PASS** | 5.5s |
| close-every-chart | atomic-undo | pass | pass | pass | pass | **PASS** | 4.2s |
| scatter-mw-logp-colored-regression | combine | pass | pass | pass | pass | **PASS** | 4.6s |
| add-then-close-sequenced | combine | pass | pass | pass | pass | **PASS** | 8.8s |
| replace-scatter-with-heatmap | combine | pass | pass | pass | pass | **PASS** | 8.6s |
| wrong-tool-trap-chem-space | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 14.9s |
| wrong-tool-trap-deprecated-shorthand | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.2s |
| wrong-tool-trap-fromtype-direct | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.0s |
| reset-to-just-grid | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.1s |
| remove-every-chart-fresh-start | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.0s |
| change-color-column | atomic-do | pass | pass | pass | pass | **PASS** | 8.7s |

## Failed prompts (deep)

### add-histogram-split-by-class — atomic-do

**Prompt:** Show me the distribution of activities per compound class
**Intent:** Distribution per category. Histogram with split column is canonical; box plot is also valid (shows distribution per group). Accept either. Reject bar chart and bare property names.

**Rubric:**

- path: PASS
- helpers: FAIL (`view\.addViewer\(\s*DG\.VIEWER\.(HISTOGRAM|BOX_PLOT)`)
- symbols: FAIL (`DG\.VIEWER\.(HISTOGRAM|BOX_PLOT)`, `(value|category)ColumnName`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`view\.addViewer\(\s*DG\.VIEWER\.(HISTOGRAM|BOX_PLOT)`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`DG\.VIEWER\.(HISTOGRAM|BOX_PLOT)`, `(value|category)ColumnName`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Inspect available columns
const cols = t.columns.toList().map(c => `${c.name} (${c.type})`);
return ui.divText(cols.join('\n'));
```

<details><summary>Full response text</summary>

```
I'll inspect the current table to find the relevant columns automatically.

```datagrok-exec
// Inspect available columns
const cols = t.columns.toList().map(c => `${c.name} (${c.type})`);
return ui.divText(cols.join('\n'));
```

Could you tell me the names of the **activity** column and the **compound class** column in your table? Once I know those, I'll add a box plot (distribution of activity per class) and optionally a histogram split by class.

Alternatively, share the column names directly and I'll build it right away! For example:
- Activity column: `activity`, `IC50`, `pIC50`, `potency`, etc.
- Compound class column: `class`, `series`, `category`, `compound_class`, etc.
```

</details>

## Suggested skill improvements

- Expected helper `view\.addViewer\(\s*DG\.VIEWER\.(HISTOGRAM|BOX_PLOT)` not reached for: add-histogram-split-by-class. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- `wrong-tool-trap-fromtype-direct`: passed rubric but emitted no exec blocks (expected_path may not have required one).
