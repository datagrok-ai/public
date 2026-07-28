# Eval report — df-and-columns

Started: 2026-05-14T09:57:21.422Z
Duration: 41.2s
Model: `sonnet`
System prompt size: 31,506 chars
Inlined skills: `datagrok-exec`, `df-and-columns`

## Summary

- **Total prompts:** 18
- **Passed:** 15
- **Failed:** 3
- **Pass rate:** 83.3%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 6 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 2 | 3 | 67% |
| reset-or-cleanup | 1 | 2 | 50% |
| wrong-tool-trap | 2 | 3 | 67% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`df-and-columns`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| find-molecule-col | atomic-do | pass | pass | pass | pass | **PASS** | 12.6s |
| describe-activity-col | atomic-do | pass | pass | pass | pass | **PASS** | 11.1s |
| add-logmw-from-mw | atomic-do | pass | pass | pass | pass | **PASS** | 7.7s |
| add-empty-float-with-units | atomic-do | pass | pass | pass | pass | **PASS** | 5.4s |
| set-friendly-name | atomic-do | pass | pass | pass | pass | **PASS** | 4.1s |
| color-coding-linear | atomic-do | pass | pass | pass | pass | **PASS** | 5.0s |
| remove-just-added | atomic-undo | pass | pass | pass | pass | **PASS** | 4.8s |
| clear-friendly-name | atomic-undo | pass | pass | pass | pass | **PASS** | 4.3s |
| rename-back | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| color-coding-off | atomic-undo | pass | pass | pass | pass | **PASS** | 5.2s |
| describe-and-add | combine | pass | pass | pass | pass | **PASS** | 25.9s |
| clone-filtered-subset | combine | pass | pass | pass | pass | **PASS** | 16.0s |
| find-and-describe-numeric | combine | pass | fail | pass | pass | **FAIL** | 11.6s |
| wrong-tool-trap-semtype | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.0s |
| wrong-tool-trap-float-literal | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 4.8s |
| wrong-tool-trap-color-tag | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 7.0s |
| reset-all-meta | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 9.5s |
| remove-session-calculated | reset-or-cleanup | pass | pass | fail | pass | **FAIL** | 9.8s |

## Failed prompts (deep)

### find-and-describe-numeric — combine

**Prompt:** summarize all numeric columns
**Intent:** Iterate t.columns.numerical or findColumns({type:'numerical'}), then describeColumn each.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.describeColumn`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.describeColumn`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const results = {};
for (const col of t.columns.numerical) {
  const s = col.stats;
  results[col.name] = {
    type: col.type,
    count: s.valueCount,
    missing: s.missingValueCount,
    min: +s.min.toFixed(4),
    max: +s.max.toFixed(4),
    avg: +s.avg.toFixed(4),
    stdev: +s.stdev.toFixed(4),
    median: +s.med.toFixed(4),
  };
}
return ui.tableFromMap(results);
```

<details><summary>Full response text</summary>

```

```datagrok-exec
const results = {};
for (const col of t.columns.numerical) {
  const s = col.stats;
  results[col.name] = {
    type: col.type,
    count: s.valueCount,
    missing: s.missingValueCount,
    min: +s.min.toFixed(4),
    max: +s.max.toFixed(4),
    avg: +s.avg.toFixed(4),
    stdev: +s.stdev.toFixed(4),
    median: +s.med.toFixed(4),
  };
}
return ui.tableFromMap(results);
```

Here's a summary of all numeric columns in your table. Each column shows:

- **count** — non-null value count
- **missing** — missing value count
- **min / max** — range
- **avg** — mean
- **stdev** — standard deviation
- **median** — middle value (Q2)
```

</details>

### wrong-tool-trap-float-literal — wrong-tool-trap

**Prompt:** make this a float column called score
**Intent:** Trap: Claude must NOT use the literal 'float'. The float type constant is the string 'double' — always use DG.COLUMN_TYPE.FLOAT.
**Trap:** DG.COLUMN_TYPE.FLOAT === 'double'. The string 'float' is not a valid column type.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.addColumn`)
- symbols: FAIL (`DG\.COLUMN_TYPE\.FLOAT`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.addColumn`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`DG\.COLUMN_TYPE\.FLOAT`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
return ui.tableFromMap(Object.fromEntries(t.columns.toList().map(c => [c.name, c.type])));
```

<details><summary>Full response text</summary>

```
I need to see the current table to understand the context. Could you clarify which column you'd like to convert to a float column named "score"? Or if you'd like me to look at the current table first:

```datagrok-exec
return ui.tableFromMap(Object.fromEntries(t.columns.toList().map(c => [c.name, c.type])));
```
```

</details>

### remove-session-calculated — reset-or-cleanup

**Prompt:** remove all calculated columns we added this session
**Intent:** Claude has to reason about which columns are custom — there is no session log. The right heuristic is: columns with a formula tag (meta.formula or tag '.formula'). Filter columns by formula presence, then removeColumns. Note this prompt requires Claude to make a judgement call — it may fail at eval time and the result will tell us whether we need a session-tracking helper.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`formula|tags|meta\.formula`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`formula|tags|meta\.formula`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Find and remove all calculated columns added this session
const calcCols = [...t.columns].filter(col => col.isCalculated);
const names = calcCols.map(c => c.name);

if (names.length === 0) {
  return ui.divText('No calculated columns found.');
}

grokky.removeColumns(t, names);
return ui.divText(`Removed ${names.length} calculated column(s): ${names.join(', ')}`);
```

<details><summary>Full response text</summary>

```
```datagrok-exec
// Find and remove all calculated columns added this session
const calcCols = [...t.columns].filter(col => col.isCalculated);
const names = calcCols.map(c => c.name);

if (names.length === 0) {
  return ui.divText('No calculated columns found.');
}

grokky.removeColumns(t, names);
return ui.divText(`Removed ${names.length} calculated column(s): ${names.join(', ')}`);
```

This will find all calculated (formula-driven) columns in the current table and remove them. I'll let you know how many were found and removed.
```

</details>

## Suggested skill improvements

- Expected helper `grokky.describeColumn` not reached for: find-and-describe-numeric. Consider tightening the skill's routing table or example density for this case.
- Expected helper `grokky.addColumn` not reached for: wrong-tool-trap-float-literal. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
