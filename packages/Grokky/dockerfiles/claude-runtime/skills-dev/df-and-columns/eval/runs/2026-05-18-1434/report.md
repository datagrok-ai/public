# Eval report — df-and-columns

Started: 2026-05-18T11:34:30.094Z
Duration: 30.1s
Model: `sonnet`
System prompt size: 36,054 chars
Inlined skills: `datagrok-exec`, `df-and-columns`

## Summary

- **Total prompts:** 18
- **Passed:** 16
- **Failed:** 2
- **Pass rate:** 88.9%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 5 | 6 | 83% |
| atomic-undo | 4 | 4 | 100% |
| combine | 2 | 3 | 67% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 30.1s wall
- Mean: 5.8s · Median: 4.4s · p95: 11.4s
- Slowest 3: `add-logmw-from-mw` (11.4s), `describe-and-add` (11.2s), `describe-activity-col` (11.0s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`df-and-columns`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| find-molecule-col | atomic-do | pass | pass | pass | pass | **PASS** | 4.5s |
| describe-activity-col | atomic-do | pass | pass | pass | pass | **PASS** | 11.0s |
| add-logmw-from-mw | atomic-do | pass | fail | fail | pass | **FAIL** | 11.4s |
| add-empty-float-with-units | atomic-do | pass | pass | pass | pass | **PASS** | 6.0s |
| set-friendly-name | atomic-do | pass | pass | pass | pass | **PASS** | 4.4s |
| color-coding-linear | atomic-do | pass | pass | pass | pass | **PASS** | 4.3s |
| remove-just-added | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| clear-friendly-name | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| rename-back | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| color-coding-off | atomic-undo | pass | pass | pass | pass | **PASS** | 3.6s |
| describe-and-add | combine | pass | fail | pass | pass | **FAIL** | 11.2s |
| clone-filtered-subset | combine | pass | pass | pass | pass | **PASS** | 4.1s |
| find-and-describe-numeric | combine | pass | pass | pass | pass | **PASS** | 6.2s |
| wrong-tool-trap-semtype | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.8s |
| wrong-tool-trap-float-literal | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.2s |
| wrong-tool-trap-color-tag | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.2s |
| reset-all-meta | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 3.8s |
| remove-session-calculated | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 7.3s |

## Failed prompts (deep)

### add-logmw-from-mw — atomic-do

**Prompt:** add a column 'logMW' computed from MW
**Intent:** Formula-driven column. Claude should reach for t.columns.addNewCalculated('logMW', '${MW}...', ...) — see datagrok-calc-column. NOT a static addNewFloat + manual loop.

**Rubric:**

- path: PASS
- helpers: FAIL (`addNewCalculated\(['"]logMW|addNewCalculated\(['"]\w+['"], *['"]\$\{MW`)
- symbols: FAIL (`Log10|Log\(`, `\$\{MW\}|\$\{mw\}`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`addNewCalculated\(['"]logMW|addNewCalculated\(['"]\w+['"], *['"]\$\{MW`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`Log10|Log\(`, `\$\{MW\}|\$\{mw\}`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
const mwCol = t.col('MW') ?? t.columns.firstWhere(c => c.name.toLowerCase().includes('mw'));
if (!mwCol) throw new Error('No MW column found in the table.');

const logMWName = t.columns.getUnusedName('logMW');
const logMWCol = t.columns.addNewFloat(logMWName);
logMWCol.init(i => {
  const v = mwCol.get(i);
  return (v != null && v > 0) ? Math.log10(v) : null;
});
return ui.divText(`Added column "${logMWName}" = log10(${mwCol.name})`);
```

<details><summary>Full response text</summary>

```
```datagrok-exec
const mwCol = t.col('MW') ?? t.columns.firstWhere(c => c.name.toLowerCase().includes('mw'));
if (!mwCol) throw new Error('No MW column found in the table.');

const logMWName = t.columns.getUnusedName('logMW');
const logMWCol = t.columns.addNewFloat(logMWName);
logMWCol.init(i => {
  const v = mwCol.get(i);
  return (v != null && v > 0) ? Math.log10(v) : null;
});
return ui.divText(`Added column "${logMWName}" = log10(${mwCol.name})`);
```

This adds a new `logMW` column (float) where each value is **log₁₀(MW)**. Rows with missing or non-positive MW values are set to `null`.
```

</details>

### describe-and-add — combine

**Prompt:** describe the molecule column and also add a Ki column with semType Ki
**Intent:** Two operations in one exec block: find molecule column via bySemType, read its stats; addNewFloat('Ki') and set col.semType = DG.SEMTYPE.Ki.

**Rubric:**

- path: PASS
- helpers: FAIL (`bySemType\(DG\.SEMTYPE\.MOLECULE\).*addNewFloat\(['"]Ki|bySemType\(DG\.SEMTYPE\.MOLECULE\).*addNew\(['"]Ki`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`bySemType\(DG\.SEMTYPE\.MOLECULE\).*addNewFloat\(['"]Ki|bySemType\(DG\.SEMTYPE\.MOLECULE\).*addNew\(['"]Ki`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
// Describe the molecule column
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!molCol) return ui.divText('No molecule column found.');

const s = molCol.stats;
const info = {
  name: molCol.name,
  type: molCol.type,
  semType: molCol.semType ?? '(none)',
  length: molCol.length,
  missing: s.missingValueCount,
  unique: s.uniqueCount,
};

// Add a Ki column with semType Ki
const kiName = t.columns.getUnusedName('Ki');
const kiCol = t.columns.addNewFloat(kiName);
kiCol.semType = DG.SEMTYPE.Ki;

return ui.tableFromMap({...info, '--- added column ---': '', 'Ki column name': kiCol.name, 'Ki semType': kiCol.semType});
```

<details><summary>Full response text</summary>

```
```datagrok-exec
// Describe the molecule column
const molCol = t.columns.bySemType(DG.SEMTYPE.MOLECULE);
if (!molCol) return ui.divText('No molecule column found.');

const s = molCol.stats;
const info = {
  name: molCol.name,
  type: molCol.type,
  semType: molCol.semType ?? '(none)',
  length: molCol.length,
  missing: s.missingValueCount,
  unique: s.uniqueCount,
};

// Add a Ki column with semType Ki
const kiName = t.columns.getUnusedName('Ki');
const kiCol = t.columns.addNewFloat(kiName);
kiCol.semType = DG.SEMTYPE.Ki;

return ui.tableFromMap({...info, '--- added column ---': '', 'Ki column name': kiCol.name, 'Ki semType': kiCol.semType});
```

Here's what the code does:

1. **Describes the molecule column** — finds it by `DG.SEMTYPE.MOLECULE` and reports:
   - Name, type, semantic type
   - Total row count, missing values, unique values

2. **Adds a `Ki` column** — creates a new `float` column with a collision-safe name and sets its semantic type to `DG.SEMTYPE.Ki`.

The Ki column is now empty (all `NaN`/missing) and ready for you to populate with values. Let me know if you'd like to fill it in or set units/format!
```

</details>

## Suggested skill improvements

- Expected helper `addNewCalculated\(['"]logMW|addNewCalculated\(['"]\w+['"], *['"]\$\{MW` not reached for: add-logmw-from-mw. Consider tightening the skill's routing table or example density for this case.
- Expected helper `bySemType\(DG\.SEMTYPE\.MOLECULE\).*addNewFloat\(['"]Ki|bySemType\(DG\.SEMTYPE\.MOLECULE\).*addNew\(['"]Ki` not reached for: describe-and-add. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
