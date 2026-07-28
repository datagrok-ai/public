# Eval report — selection

Started: 2026-05-18T11:35:11.042Z
Duration: 30.0s
Model: `sonnet`
System prompt size: 34,415 chars
Inlined skills: `datagrok-exec`, `selection`

## Summary

- **Total prompts:** 18
- **Passed:** 17
- **Failed:** 1
- **Pass rate:** 94.4%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 5 | 5 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 4 | 75% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 30.0s wall
- Mean: 5.7s · Median: 4.5s · p95: 11.2s
- Slowest 3: `reset-select-all-molecule-column` (11.2s), `describe-current-selection` (10.7s), `select-top-5-by-activity` (10.2s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`selection`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| select-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 5.7s |
| select-rows-0-10 | atomic-do | pass | pass | pass | pass | **PASS** | 4.9s |
| select-top-5-by-activity | atomic-do | pass | pass | pass | pass | **PASS** | 10.2s |
| select-ic50-under-100 | atomic-do | pass | pass | pass | pass | **PASS** | 7.3s |
| select-all-rows | atomic-do | pass | pass | pass | pass | **PASS** | 9.2s |
| clear-the-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| deselect-everything | atomic-undo | pass | pass | pass | pass | **PASS** | 4.3s |
| invert-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| remove-rows-5-and-7 | atomic-undo | pass | pass | pass | pass | **PASS** | 4.8s |
| show-selected-as-new-table | combine | pass | pass | pass | pass | **PASS** | 4.5s |
| select-all-visible-rows | combine | pass | pass | pass | pass | **PASS** | 4.1s |
| filter-to-selected-rows | combine | pass | pass | pass | pass | **PASS** | 4.3s |
| wrong-tool-trap-clear-via-setall-true | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.8s |
| wrong-tool-trap-select-all-not-clear | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.3s |
| wrong-tool-trap-highlight-row-5-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.8s |
| reset-deselect-and-default | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.0s |
| reset-select-all-molecule-column | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 11.2s |
| describe-current-selection | combine | pass | fail | pass | pass | **FAIL** | 10.7s |

## Failed prompts (deep)

### describe-current-selection — combine

**Prompt:** describe what's currently selected
**Intent:** Hand-roll a report: count, capped list of indices, optional sample row from the first selected index. Must NOT read t.filter — wrong mask.

**Rubric:**

- path: PASS
- helpers: FAIL (`selection\.trueCount|selection\.getSelectedIndexes`)
- symbols: PASS
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`selection\.trueCount|selection\.getSelectedIndexes`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper.

**Emitted exec blocks:**

```javascript
const sel = t.selection;
const all = sel.getSelectedIndexes();
const MAX_INDEXES = 50;
const MAX_SAMPLE_COLS = 12;
const indexes = Array.from(all.slice(0, MAX_INDEXES));

const out = {
  count: sel.trueCount,
  total: t.rowCount,
  indexes,
  currentRowIdx: t.currentRowIdx,
  sample: undefined,
};

if (all.length > 0) {
  const firstIdx = all[0];
  const cols = t.columns;
  const sample = {};
  const colCap = Math.min(cols.length, MAX_SAMPLE_COLS);
  for (let c = 0; c < colCap; c++) {
    const col = cols.byIndex(c);
    sample[col.name] = col.get(firstIdx);
  }
  out.sample = sample;
}
return ui.tableFromMap(out);
```

<details><summary>Full response text</summary>

```
```datagrok-exec
const sel = t.selection;
const all = sel.getSelectedIndexes();
const MAX_INDEXES = 50;
const MAX_SAMPLE_COLS = 12;
const indexes = Array.from(all.slice(0, MAX_INDEXES));

const out = {
  count: sel.trueCount,
  total: t.rowCount,
  indexes,
  currentRowIdx: t.currentRowIdx,
  sample: undefined,
};

if (all.length > 0) {
  const firstIdx = all[0];
  const cols = t.columns;
  const sample = {};
  const colCap = Math.min(cols.length, MAX_SAMPLE_COLS);
  for (let c = 0; c < colCap; c++) {
    const col = cols.byIndex(c);
    sample[col.name] = col.get(firstIdx);
  }
  out.sample = sample;
}
return ui.tableFromMap(out);
```

This block reports:
- **count** — how many rows are currently selected
- **total** — total rows in the table
- **indexes** — the row indices that are selected (up to the first 50)
- **currentRowIdx** — the focused/current row
- **sample** — column values from the first selected row (up to 12 columns)
```

</details>

## Suggested skill improvements

- Expected helper `selection\.trueCount|selection\.getSelectedIndexes` not reached for: describe-current-selection. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
