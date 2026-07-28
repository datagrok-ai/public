# Eval report — df-and-columns

Started: 2026-05-18T11:40:20.163Z
Duration: 34.6s
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

- Total: 34.6s wall
- Mean: 6.6s · Median: 5.1s · p95: 14.5s
- Slowest 3: `add-logmw-from-mw` (14.5s), `describe-activity-col` (13.0s), `describe-and-add` (11.5s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`df-and-columns`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| find-molecule-col | atomic-do | pass | pass | pass | pass | **PASS** | 5.3s |
| describe-activity-col | atomic-do | pass | pass | pass | pass | **PASS** | 13.0s |
| add-logmw-from-mw | atomic-do | pass | fail | pass | pass | **FAIL** | 14.5s |
| add-empty-float-with-units | atomic-do | pass | pass | pass | pass | **PASS** | 9.2s |
| set-friendly-name | atomic-do | pass | pass | pass | pass | **PASS** | 5.1s |
| color-coding-linear | atomic-do | pass | pass | pass | pass | **PASS** | 4.6s |
| remove-just-added | atomic-undo | pass | pass | pass | pass | **PASS** | 4.6s |
| clear-friendly-name | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| rename-back | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| color-coding-off | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| describe-and-add | combine | pass | pass | fail | pass | **FAIL** | 11.5s |
| clone-filtered-subset | combine | pass | pass | pass | pass | **PASS** | 3.9s |
| find-and-describe-numeric | combine | pass | pass | pass | pass | **PASS** | 6.7s |
| wrong-tool-trap-semtype | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 6.0s |
| wrong-tool-trap-float-literal | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.7s |
| wrong-tool-trap-color-tag | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.1s |
| reset-all-meta | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.9s |
| remove-session-calculated | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 10.0s |

## Failed prompts (deep)

### add-logmw-from-mw — atomic-do

**Prompt:** add a column 'logMW' computed from MW
**Intent:** Calculated column. Preferred path is the formula DSL via addNewCalculated('logMW', '${MW}...') (datagrok-calc-column owns formulas). In isolated eval the static fallback (addNewFloat('logMW') + col.init(i => Math.log10(...))) is also accepted. Reject the per-row imperative for-loop.

**Rubric:**

- path: PASS
- helpers: FAIL (`addNewCalculated\(['"]logMW|addNewCalculated\(['"]\w+['"], *['"]\$\{MW|addNewFloat\(['"]logMW[\s\S]*\.init\([\s\S]*Math\.log10`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`addNewCalculated\(['"]logMW|addNewCalculated\(['"]\w+['"], *['"]\$\{MW|addNewFloat\(['"]logMW[\s\S]*\.init\([\s\S]*Math\.log10`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const mwCol = t.col('MW') ?? t.columns.firstWhere(c => c.name.toLowerCase() === 'mw');
if (!mwCol) throw new Error('Column "MW" not found in the table.');

const logMW = t.columns.addNewFloat(t.columns.getUnusedName('logMW'));
logMW.init(i => {
  const v = mwCol.get(i);
  return (v == null || v <= 0) ? null : Math.log10(v);
});
logMW.meta.description = 'log10(MW)';
return ui.divText(`Added column "${logMW.name}" (log10 of "${mwCol.name}")`);
```

<details><summary>Full response text</summary>

```
```datagrok-exec
const mwCol = t.col('MW') ?? t.columns.firstWhere(c => c.name.toLowerCase() === 'mw');
if (!mwCol) throw new Error('Column "MW" not found in the table.');

const logMW = t.columns.addNewFloat(t.columns.getUnusedName('logMW'));
logMW.init(i => {
  const v = mwCol.get(i);
  return (v == null || v <= 0) ? null : Math.log10(v);
});
logMW.meta.description = 'log10(MW)';
return ui.divText(`Added column "${logMW.name}" (log10 of "${mwCol.name}")`);
```

This adds a new `logMW` column computed as **log₁₀(MW)**:
- Rows where MW is missing or ≤ 0 are set to `null`.
- The column type is `float`.
- The name is collision-safe (e.g. `logMW (1)` if `logMW` already exists).
```

</details>

### describe-and-add — combine

**Prompt:** describe the molecule column and also add a Ki column with semType Ki
**Intent:** Two operations in one exec block: find molecule column via bySemType, read its stats; addNewFloat('Ki') and set col.semType = DG.SEMTYPE.Ki.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`addNewFloat\(['"]Ki|addNew\(['"]Ki`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`addNewFloat\(['"]Ki|addNew\(['"]Ki`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

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

// Add Ki column with semType Ki
const kiName = t.columns.getUnusedName('Ki');
const kiCol = t.columns.addNewFloat(kiName);
kiCol.semType = DG.SEMTYPE.Ki;
kiCol.meta.units = 'nM';

return ui.divV([
  ui.h2('Molecule Column Description'),
  ui.tableFromMap(info),
  ui.divText(`✓ Added column "${kiCol.name}" (type: ${kiCol.type}, semType: ${kiCol.semType}, units: ${kiCol.meta.units})`),
]);
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

// Add Ki column with semType Ki
const kiName = t.columns.getUnusedName('Ki');
const kiCol = t.columns.addNewFloat(kiName);
kiCol.semType = DG.SEMTYPE.Ki;
kiCol.meta.units = 'nM';

return ui.divV([
  ui.h2('Molecule Column Description'),
  ui.tableFromMap(info),
  ui.divText(`✓ Added column "${kiCol.name}" (type: ${kiCol.type}, semType: ${kiCol.semType}, units: ${kiCol.meta.units})`),
]);
```

Here's what was done:

### 🔬 Molecule Column
The block finds the first column with `semType = Molecule` and reports:
- **Name, type, semType** — to confirm it's the right column
- **Length** — total row count
- **Missing values** — null/empty count
- **Unique values** — distinct structure count

### ➕ Ki Column Added
A new **`Ki`** column was added with:
- **Type:** `float` (`DG.COLUMN_TYPE.FLOAT`)
- **SemType:** `DG.SEMTYPE.Ki`
- **Units:** `nM`

The name is collision-safe — if `Ki` already existed, it would have been named `Ki (1)`, etc.
```

</details>

## Suggested skill improvements

- Expected helper `addNewCalculated\(['"]logMW|addNewCalculated\(['"]\w+['"], *['"]\$\{MW|addNewFloat\(['"]logMW[\s\S]*\.init\([\s\S]*Math\.log10` not reached for: add-logmw-from-mw. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
