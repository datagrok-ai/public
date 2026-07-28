# Eval report — filtering

Started: 2026-05-14T10:21:59.817Z
Duration: 40.4s
Model: `sonnet`
System prompt size: 34,406 chars
Inlined skills: `datagrok-exec`, `filtering`

## Summary

- **Total prompts:** 18
- **Passed:** 14
- **Failed:** 4
- **Pass rate:** 77.8%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 6 | 100% |
| atomic-undo | 3 | 4 | 75% |
| combine | 2 | 3 | 67% |
| reset-or-cleanup | 1 | 2 | 50% |
| wrong-tool-trap | 2 | 3 | 67% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`filtering`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| filter-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 6.5s |
| filter-mw-lt-500 | atomic-do | pass | pass | pass | pass | **PASS** | 6.8s |
| filter-substructure-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 5.4s |
| filter-categorical-in-set | atomic-do | pass | pass | pass | pass | **PASS** | 4.7s |
| filter-name-contains | atomic-do | pass | pass | pass | pass | **PASS** | 3.7s |
| filter-smarts-pattern | atomic-do | pass | pass | pass | pass | **PASS** | 4.5s |
| clear-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| show-all-rows | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| remove-substructure-only | atomic-undo | pass | fail | fail | pass | **FAIL** | 9.1s |
| invert-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 3.6s |
| filter-mw-and-logp | combine | pass | pass | pass | pass | **PASS** | 5.8s |
| filtered-as-new-table | combine | pass | pass | pass | pass | **PASS** | 5.7s |
| cross-skill-selected-then-describe | combine | pass | fail | pass | pass | **FAIL** | 26.9s |
| wrong-tool-trap-remove-as-delete | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.4s |
| wrong-tool-trap-smiles-as-molblock | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.4s |
| wrong-tool-trap-hide-vs-drop | wrong-tool-trap | fail | fail | fail | pass | **FAIL** | 9.0s |
| reset-clear-all-filters-and-selections | reset-or-cleanup | pass | pass | pass | fail | **FAIL** | 4.3s |
| go-back-to-seeing-every-row | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.3s |

## Failed prompts (deep)

### remove-substructure-only — atomic-undo

**Prompt:** remove the substructure filter only
**Intent:** Only the substructure filter should go away — other filters stay. Acceptable: filterSubstructure(view, col, '') (empty query → clear path), or fg.remove() targeting the substructure filter widget. NOT clearFilter(view) which disables everything.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.filterSubstructure`)
- symbols: FAIL (`filterSubstructure.*['"]['"]|molBlock: ['"]['"]|getStates`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.filterSubstructure`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`filterSubstructure.*['"]['"]|molBlock: ['"]['"]|getStates`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const fg = view.getFiltersGroup({createDefaultFilters: false});
const filters = fg.filters;
for (const f of filters) {
  if (f.filterType === DG.FILTER_TYPE.SUBSTRUCTURE) {
    fg.remove(f);
    break;
  }
}
```

<details><summary>Full response text</summary>

```
```datagrok-exec
const fg = view.getFiltersGroup({createDefaultFilters: false});
const filters = fg.filters;
for (const f of filters) {
  if (f.filterType === DG.FILTER_TYPE.SUBSTRUCTURE) {
    fg.remove(f);
    break;
  }
}
```

This removes the substructure filter from the filter panel while leaving all other active filters (range, categorical, free-text, etc.) intact.
```

</details>

### cross-skill-selected-then-describe — combine

**Prompt:** filter to selected molecules then describe them
**Intent:** Cross-skill: needs to bridge selection (other skill) and filter (this skill). Acceptable path: filterByPredicate(df, i => df.selection.get(i)) or filteredDf using a clone with the selection mask. Pass if Claude reaches across cleanly.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.filterByPredicate|grokky.filteredDf`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.filterByPredicate|grokky.filteredDf`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
// Step 1: Copy selection into the filter
const sel = t.selection;
if (sel.trueCount === 0) {
  return ui.divText('⚠️ No molecules are selected. Please select some rows first (e.g. via lasso or Shift+click), then run this again.');
}

// Apply the selection mask as the filter
t.filter.copyFrom(sel);

// Step 2: Find the molecule column
const molCol = t.columns.toList().find(c =>
  c.semType === DG.SEMTYPE.MOLECULE || c.semType === 'Molecule'
);

// Step 3: Build a summary
const n = t.filter.trueCount;
const total = t.rowCount;

const lines = [`**${n} of ${total} rows** now visible (filtered to your selection).`];

// Numeric column stats on filtered rows
const numCols = t.columns.toList().filter(c => c.type === DG.TYPE.FLOAT || c.type === DG.TYPE.INT);
const indices = [...t.filter.getSelectedIndexes()];

const statRows = numCols.slice(0, 8).map(col => {
  const vals = indices.map(i => col.get(i)).filter(v => v !== null && !isNaN(v));
  if (!vals.length) return null;
  const min = Math.min(...vals).toFixed(2);
  const max = Math.max(...vals).toFixed(2);
  const mean = (vals.reduce((a, b) => a + b, 0) / vals.length).toFixed(2);
  return { Column: col.name, Min: min, Mean: mean, Max: max, Count: vals.length };
}).filter(Boolean);

const catCols = t.columns.toList().filter(c => c.type === DG.TYPE.STRING && c !== molCol);
const catSummary = catCols.slice(0, 3).map(col => {
  const freq = {};
  indices.forEach(i => { const v = col.get(i) ?? '(null)'; freq[v] = (freq[v] || 0) + 1; });
  const top = Object.entries(freq).sort((a,b) => b[1]-a[1]).slice(0, 3).map(([k,v]) => `${k} (${v})`).join(', ');
  return `**${col.name}:** ${top}`;
});

const container = ui.divV([
  ui.h2(`Selection Summary — ${n} molecule${n !== 1 ? 's' : ''}`),
  molCol ? ui.divText(`🔬 Molecule column: ${molCol.name}`) : ui.divText('⚠️ No molecule column detected'),
  ui.divText(''),
]);

if (statRows.length) {
  const tbl = ui.tableFromMap(Object.fromEntries(statRows.map(r =>
    [r.Column, `min ${r.Min}  |  mean ${r.Mean}  |  max ${r.Max}`]
  )));
  container.append(ui.h3('Numeric columns'));
  container.append(tbl);
}

if (catSummary.length) {
  container.append(ui.h3('Categorical columns (top values)'));
  catSummary.forEach(s => container.append(ui.markdown(s)));
}

return container;
```

<details><summary>Full response text</summary>

```
Here's the plan:
1. Filter the table to show only the currently **selected** rows
2. Describe the visible molecules (summary statistics or structure info)

Let me first filter to selected rows, then describe them:

```datagrok-exec
// Step 1: Copy selection into the filter
const sel = t.selection;
if (sel.trueCount === 0) {
  return ui.divText('⚠️ No molecules are selected. Please select some rows first (e.g. via lasso or Shift+click), then run this again.');
}

// Apply the selection mask as the filter
t.filter.copyFrom(sel);

// Step 2: Find the molecule column
const molCol = t.columns.toList().find(c =>
  c.semType === DG.SEMTYPE.MOLECULE || c.semType === 'Molecule'
);

// Step 3: Build a summary
const n = t.filter.trueCount;
const total = t.rowCount;

const lines = [`**${n} of ${total} rows** now visible (filtered to your selection).`];

// Numeric column stats on filtered rows
const numCols = t.columns.toList().filter(c => c.type === DG.TYPE.FLOAT || c.type === DG.TYPE.INT);
const indices = [...t.filter.getSelectedIndexes()];

const statRows = numCols.slice(0, 8).map(col => {
  const vals = indices.map(i => col.get(i)).filter(v => v !== null && !isNaN(v));
  if (!vals.length) return null;
  const min = Math.min(...vals).toFixed(2);
  const max = Math.max(...vals).toFixed(2);
  const mean = (vals.reduce((a, b) => a + b, 0) / vals.length).toFixed(2);
  return { Column: col.name, Min: min, Mean: mean, Max: max, Count: vals.length };
}).filter(Boolean);

const catCols = t.columns.toList().filter(c => c.type === DG.TYPE.STRING && c !== molCol);
const catSummary = catCols.slice(0, 3).map(col => {
  const freq = {};
  indices.forEach(i => { const v = col.get(i) ?? '(null)'; freq[v] = (freq[v] || 0) + 1; });
  const top = Object.entries(freq).sort((a,b) => b[1]-a[1]).slice(0, 3).map(([k,v]) => `${k} (${v})`).join(', ');
  return `**${col.name}:** ${top}`;
});

const container = ui.divV([
  ui.h2(`Selection Summary — ${n} molecule${n !== 1 ? 's' : ''}`),
  molCol ? ui.divText(`🔬 Molecule column: ${molCol.name}`) : ui.divText('⚠️ No molecule column detected'),
  ui.divText(''),
]);

if (statRows.length) {
  const tbl = ui.tableFromMap(Object.fromEntries(statRows.map(r =>
    [r.Column, `min ${r.Min}  |  mean ${r.Mean}  |  max ${r.Max}`]
  )));
  container.append(ui.h3('Numeric columns'));
  container.append(tbl);
}

if (catSummary.length) {
  container.append(ui.h3('Categorical columns (top values)'));
  catSummary.forEach(s => container.append(ui.markdown(s)));
}

return container;
```

This will:
- **Filter** the table to show only your selected molecules
- **Report** numeric stats (min/mean/max) for each numeric column across the selection
- **Summarize** top categories for string columns

If nothing is selected yet, it will tell you — just select some rows first (lasso in a scatter plot, or Shift+click rows) and run again!The code above is already correct and ready to run — the `df.selection` BitSet and `DG.SEMTYPE.MOLECULE` patterns are standard Datagrok JS API. The block above will execute as-is. Just make sure you have some rows selected in your table view first!
```

</details>

### wrong-tool-trap-hide-vs-drop — wrong-tool-trap

**Prompt:** hide rows where activity is high
**Intent:** Trap: 'hide' is non-destructive — must use filter, not drop. Also polarity: 'hide where high' means filter to KEEP rows where activity is NOT high. Pass if Claude either filters non-destructively with inverted polarity, or asks for clarification on the threshold for 'high'.
**Trap:** 'hide' → non-destructive (filter), not removeWhere. Polarity must invert: hide_where_X is equivalent to keep_where_not_X.

**Rubric:**

- path: FAIL (`>=1 datagrok-exec block`)
- helpers: FAIL (`grokky.filterByPredicate|grokky.filterRows`)
- symbols: FAIL (`activity`)
- forbidden: PASS

**Hypothesis:** Claude did not emit a `datagrok-exec` block at all — likely treating the request as conversational. Skill may need a stronger "always emit a block" directive for this intent class. Claude bypassed the expected helper(s) (`grokky.filterByPredicate|grokky.filterRows`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`activity`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:** (none)

<details><summary>Full response text</summary>

```
I need a couple of details to apply the filter correctly:

1. **Which column** holds the activity values? (e.g. `activity`, `IC50`, `% inhibition`, …)
2. **What threshold** defines "high"? (e.g. > 7, > 100, top 10%, …)

Please share those and I'll filter the rows right away!
```

</details>

### reset-clear-all-filters-and-selections — reset-or-cleanup

**Prompt:** clear all filters and selections
**Intent:** Cross-skill: filter clearing is this skill; selection clearing is the next skill. Pass if Claude clears the filter correctly. Bonus if Claude also handles selection (df.selection.setAll(false)) — note that selection's clear-polarity is OPPOSITE the filter's.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: PASS
- forbidden: FAIL (`setAll\(false\)`)

**Hypothesis:** Forbidden pattern(s) matched (`setAll\(false\)`). Skill's anti-pattern callout isn't strong enough — consider making the wrong path throw, or surface a more obvious "DO NOT" callout near the relevant example.

**Emitted exec blocks:**

```javascript
grokky.clearFilter(view);
t.selection.setAll(false);
```

<details><summary>Full response text</summary>

```

```datagrok-exec
grokky.clearFilter(view);
t.selection.setAll(false);
```

Done — all filters cleared and all row selections deselected.
```

</details>

## Suggested skill improvements

- Expected helper `grokky.filterSubstructure` not reached for: remove-substructure-only. Consider tightening the skill's routing table or example density for this case.
- Expected helper `grokky.filterByPredicate|grokky.filteredDf` not reached for: cross-skill-selected-then-describe. Consider tightening the skill's routing table or example density for this case.
- Expected helper `grokky.filterByPredicate|grokky.filterRows` not reached for: wrong-tool-trap-hide-vs-drop. Consider tightening the skill's routing table or example density for this case.
- Forbidden pattern `setAll\(false\)` slipped through for: reset-clear-all-filters-and-selections. Anti-pattern callout needs to be more prominent — consider a "WRONG / RIGHT" pair right before the relevant example.

## Anomalies

- None.
