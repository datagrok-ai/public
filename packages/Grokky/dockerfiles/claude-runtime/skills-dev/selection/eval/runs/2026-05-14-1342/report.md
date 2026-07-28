# Eval report — selection

Started: 2026-05-14T10:42:12.334Z
Duration: 53.5s
Model: `sonnet`
System prompt size: 33,974 chars
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
| combine | 4 | 4 | 100% |
| reset-or-cleanup | 1 | 2 | 50% |
| wrong-tool-trap | 3 | 3 | 100% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`selection`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| select-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 6.4s |
| select-rows-0-10 | atomic-do | pass | pass | pass | pass | **PASS** | 4.7s |
| select-top-5-by-activity | atomic-do | pass | pass | pass | pass | **PASS** | 17.1s |
| select-ic50-under-100 | atomic-do | pass | pass | pass | pass | **PASS** | 10.0s |
| select-all-rows | atomic-do | pass | pass | pass | pass | **PASS** | 8.8s |
| clear-the-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 4.1s |
| deselect-everything | atomic-undo | pass | pass | pass | pass | **PASS** | 4.7s |
| invert-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 4.6s |
| remove-rows-5-and-7 | atomic-undo | pass | pass | pass | pass | **PASS** | 3.9s |
| show-selected-as-new-table | combine | pass | pass | pass | pass | **PASS** | 4.1s |
| select-all-visible-rows | combine | pass | pass | pass | pass | **PASS** | 5.6s |
| filter-to-selected-rows | combine | pass | pass | pass | pass | **PASS** | 4.5s |
| wrong-tool-trap-clear-via-setall-true | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.6s |
| wrong-tool-trap-select-all-not-clear | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.9s |
| wrong-tool-trap-highlight-row-5-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.1s |
| reset-deselect-and-default | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.1s |
| reset-select-all-molecule-column | reset-or-cleanup | pass | pass | fail | pass | **FAIL** | 11.3s |
| describe-current-selection | combine | pass | pass | pass | pass | **PASS** | 30.8s |

## Failed prompts (deep)

### reset-select-all-molecule-column — reset-or-cleanup

**Prompt:** select all the molecule column
**Intent:** TEACHING PROMPT: columns don't have selection state — selection is per-row. The right answer is to select all rows (because every row has a molecule cell). Pass if Claude selects all rows. Bonus if Claude clarifies that selection is a row concept.
**Trap:** Selection is row-shaped, not column-shaped. 'Select all the X column' must map to 'select all rows' — there is no column selection BitSet.

**Rubric:**

- path: PASS
- helpers: PASS
- symbols: FAIL (`selectAll|setAll\(true\)`)
- forbidden: PASS

**Hypothesis:** Expected symbol(s) absent (`selectAll|setAll\(true\)`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Find the molecule column (semantic type 'Molecule') and select all rows that have a value in it
const molCol = Array.from(t.columns).find(c => c.semType === 'Molecule' || c.semType === DG.SEMTYPE.MOLECULE);

if (!molCol) {
  return ui.divText('No molecule column found in the current table.');
}

grokky.selectRows(t, (i) => molCol.get(i) != null && molCol.get(i) !== '');

return ui.tableFromMap({
  'Molecule column': molCol.name,
  'Rows selected': t.selection.trueCount,
  'Total rows': t.rowCount
});
```

<details><summary>Full response text</summary>

```
I'll select all rows where the molecule column has a value. Let me first identify the molecule column in the current table.

```datagrok-exec
// Find the molecule column (semantic type 'Molecule') and select all rows that have a value in it
const molCol = Array.from(t.columns).find(c => c.semType === 'Molecule' || c.semType === DG.SEMTYPE.MOLECULE);

if (!molCol) {
  return ui.divText('No molecule column found in the current table.');
}

grokky.selectRows(t, (i) => molCol.get(i) != null && molCol.get(i) !== '');

return ui.tableFromMap({
  'Molecule column': molCol.name,
  'Rows selected': t.selection.trueCount,
  'Total rows': t.rowCount
});
```

This selects all rows in the molecule column (any row with a non-null, non-empty molecule value). The result summary shows which column was found and how many rows were selected.
```

</details>

## Suggested skill improvements

- No structural issues detected. Either the eval passed broadly, or the failures are one-offs requiring case-by-case review.

## Anomalies

- None.
