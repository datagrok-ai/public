# Eval report — selection

Started: 2026-05-14T10:41:05.524Z
Duration: 27.7s
Model: `sonnet`
System prompt size: 33,316 chars
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

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`selection`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| select-age-gt-40 | atomic-do | pass | pass | pass | pass | **PASS** | 5.6s |
| select-rows-0-10 | atomic-do | pass | pass | pass | pass | **PASS** | 5.7s |
| select-top-5-by-activity | atomic-do | pass | pass | pass | pass | **PASS** | 8.4s |
| select-ic50-under-100 | atomic-do | pass | pass | pass | pass | **PASS** | 5.3s |
| select-all-rows | atomic-do | pass | pass | pass | pass | **PASS** | 8.7s |
| clear-the-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 4.5s |
| deselect-everything | atomic-undo | pass | pass | pass | pass | **PASS** | 4.7s |
| invert-selection | atomic-undo | pass | pass | pass | pass | **PASS** | 4.0s |
| remove-rows-5-and-7 | atomic-undo | pass | pass | pass | pass | **PASS** | 4.3s |
| show-selected-as-new-table | combine | pass | fail | fail | fail | **FAIL** | 3.6s |
| select-all-visible-rows | combine | pass | pass | pass | pass | **PASS** | 5.4s |
| filter-to-selected-rows | combine | pass | pass | pass | pass | **PASS** | 4.2s |
| wrong-tool-trap-clear-via-setall-true | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.2s |
| wrong-tool-trap-select-all-not-clear | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.8s |
| wrong-tool-trap-highlight-row-5-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 4.5s |
| reset-deselect-and-default | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 4.2s |
| reset-select-all-molecule-column | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 9.4s |
| describe-current-selection | combine | pass | pass | pass | pass | **PASS** | 3.9s |

## Failed prompts (deep)

### show-selected-as-new-table — combine

**Prompt:** show only the selected rows as a new table
**Intent:** Non-destructive new DF from selection. selectedDf(df) or df.clone(df.selection). Should NOT mutate the source via filterFromSelection (that hides rows but doesn't make a new table).

**Rubric:**

- path: PASS
- helpers: FAIL (`grokky.selectedDf`)
- symbols: FAIL (`selectedDf|clone\(.*selection`)
- forbidden: FAIL (`df\.filter\.copyFrom|filterFromSelection`)

**Hypothesis:** Claude bypassed the expected helper(s) (`grokky.selectedDf`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`selectedDf|clone\(.*selection`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill. Forbidden pattern(s) matched (`df\.filter\.copyFrom|filterFromSelection`). Skill's anti-pattern callout isn't strong enough — consider making the wrong path throw, or surface a more obvious "DO NOT" callout near the relevant example.

**Emitted exec blocks:**

```javascript
grokky.filterFromSelection(t);
```

<details><summary>Full response text</summary>

```

```datagrok-exec
grokky.filterFromSelection(t);
```

Done! The table now shows only the selected rows — the filter has been set to match your selection.
```

</details>

## Suggested skill improvements

- Expected helper `grokky.selectedDf` not reached for: show-selected-as-new-table. Consider tightening the skill's routing table or example density for this case.
- Forbidden pattern `df\.filter\.copyFrom|filterFromSelection` slipped through for: show-selected-as-new-table. Anti-pattern callout needs to be more prominent — consider a "WRONG / RIGHT" pair right before the relevant example.

## Anomalies

- None.
