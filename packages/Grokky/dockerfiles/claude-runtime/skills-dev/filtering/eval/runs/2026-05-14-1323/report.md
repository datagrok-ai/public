# Eval report — filtering

Started: 2026-05-14T10:23:41.480Z
Duration: 37.0s
Model: `sonnet`
System prompt size: 34,406 chars
Inlined skills: `datagrok-exec`, `filtering`

## Summary

- **Total prompts:** 18
- **Passed:** 17
- **Failed:** 1
- **Pass rate:** 94.4%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 6 | 6 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 2 | 2 | 100% |
| wrong-tool-trap | 2 | 3 | 67% |

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`filtering`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| filter-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 5.6s |
| filter-mw-lt-500 | atomic-do | pass | pass | pass | pass | **PASS** | 9.6s |
| filter-substructure-benzene | atomic-do | pass | pass | pass | pass | **PASS** | 5.5s |
| filter-categorical-in-set | atomic-do | pass | pass | pass | pass | **PASS** | 5.6s |
| filter-name-contains | atomic-do | pass | pass | pass | pass | **PASS** | 5.1s |
| filter-smarts-pattern | atomic-do | pass | pass | pass | pass | **PASS** | 4.8s |
| clear-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| show-all-rows | atomic-undo | pass | pass | pass | pass | **PASS** | 3.6s |
| remove-substructure-only | atomic-undo | pass | pass | pass | pass | **PASS** | 9.6s |
| invert-filter | atomic-undo | pass | pass | pass | pass | **PASS** | 4.3s |
| filter-mw-and-logp | combine | pass | pass | pass | pass | **PASS** | 5.8s |
| filtered-as-new-table | combine | pass | pass | pass | pass | **PASS** | 4.6s |
| cross-skill-selected-then-describe | combine | pass | pass | pass | pass | **PASS** | 22.3s |
| wrong-tool-trap-remove-as-delete | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.4s |
| wrong-tool-trap-smiles-as-molblock | wrong-tool-trap | pass | fail | fail | pass | **FAIL** | 9.0s |
| wrong-tool-trap-hide-vs-drop | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 5.8s |
| reset-clear-all-filters-and-selections | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 5.5s |
| go-back-to-seeing-every-row | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.5s |

## Failed prompts (deep)

### wrong-tool-trap-smiles-as-molblock — wrong-tool-trap

**Prompt:** filter by substructure CCO
**Intent:** Trap: Claude must NOT pass 'CCO' (SMILES) as molBlock in a raw substructure filter state. The fixed helper auto-detects and converts.
**Trap:** Chem:substructureFilter reads state.molBlock; raw SMILES silently matches zero rows. Auto-detect via DG.chem.isMolBlock + DG.chem.convert.

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.filterSubstructure|grokky.filterRows`)
- symbols: FAIL (`CCO`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.filterSubstructure|grokky.filterRows`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`CCO`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
// Find molecule columns in the current table
const cols = [];
for (const col of t.columns) {
  cols.push(`${col.name} [${col.type}] semType=${col.semType}`);
}
return ui.divV(cols.map(c => ui.divText(c)));
```

<details><summary>Full response text</summary>

```
Let me check the current table's columns to find the molecule column.

```datagrok-exec
// Find molecule columns in the current table
const cols = [];
for (const col of t.columns) {
  cols.push(`${col.name} [${col.type}] semType=${col.semType}`);
}
return ui.divV(cols.map(c => ui.divText(c)));
```
```

</details>

## Suggested skill improvements

- Expected helper `grokky.filterSubstructure|grokky.filterRows` not reached for: wrong-tool-trap-smiles-as-molblock. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
