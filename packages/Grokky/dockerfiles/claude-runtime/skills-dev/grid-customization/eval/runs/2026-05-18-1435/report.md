# Eval report — grid-customization

Started: 2026-05-18T11:35:18.669Z
Duration: 31.3s
Model: `sonnet`
System prompt size: 25,838 chars
Inlined skills: `datagrok-exec`, `grid-customization`

## Summary

- **Total prompts:** 19
- **Passed:** 18
- **Failed:** 1
- **Pass rate:** 94.7%

### By category

| Category | Passed | Total | Rate |
|----------|--------|-------|------|
| atomic-do | 7 | 7 | 100% |
| atomic-undo | 4 | 4 | 100% |
| combine | 3 | 3 | 100% |
| reset-or-cleanup | 1 | 2 | 50% |
| wrong-tool-trap | 3 | 3 | 100% |

### Latency

- Total: 31.3s wall
- Mean: 6.0s · Median: 4.8s · p95: 16.1s
- Slowest 3: `widen-molecule-column` (16.1s), `sort-activity-desc` (9.8s), `trap-sortbycolumns` (8.8s)

## Configuration notes

- **Test mode:** direct `@anthropic-ai/claude-agent-sdk` `query()` call, not the live container. The container is a thin Hono/WebSocket wrapper around the same SDK; we mirror `buildOptions()` from `claude-runtime/src/server.ts` so Claude sees the same `systemPrompt`/`model`/`effort`/`thinking` shape.
- **Skills inlined into system prompt:** `datagrok-exec` (the fenced-block contract — always present in prod) plus the skill under test (`grid-customization`). The four other production inlined skills (`datagrok-viewer`, `datagrok-calc-column`, `datagrok-table-ops`, `datagrok-projects`) are intentionally excluded to isolate the skill under test.
- **Tools allowed:** none. We test what Claude *writes*, not what it does with Bash/Grep. Tool use would not change the rubric outcome (rubrics check the textual exec blocks) but would balloon cost.
- **Sampling:** `model=sonnet`, `effort=low`, `thinking=disabled` — exact match to production.

## Per-prompt results

| id | category | path | helpers | symbols | forbidden | overall | latency |
|----|----------|------|---------|---------|-----------|---------|---------|
| sort-activity-desc | atomic-do | pass | pass | pass | pass | **PASS** | 9.8s |
| hide-index-column | atomic-do | pass | pass | pass | pass | **PASS** | 4.6s |
| pin-smiles-left | atomic-do | pass | pass | pass | pass | **PASS** | 4.7s |
| widen-molecule-column | atomic-do | pass | pass | pass | pass | **PASS** | 16.1s |
| color-activity-red-green | atomic-do | pass | pass | pass | pass | **PASS** | 4.8s |
| format-ic50-two-decimals | atomic-do | pass | pass | pass | pass | **PASS** | 3.8s |
| show-only-five-columns | atomic-do | pass | pass | pass | pass | **PASS** | 5.3s |
| clear-sort | atomic-undo | pass | pass | pass | pass | **PASS** | 3.5s |
| show-all-hidden | atomic-undo | pass | pass | pass | pass | **PASS** | 3.7s |
| unpin-smiles | atomic-undo | pass | pass | pass | pass | **PASS** | 3.7s |
| color-off-activity | atomic-undo | pass | pass | pass | pass | **PASS** | 4.2s |
| hit-list-multi-step | combine | pass | pass | pass | pass | **PASS** | 6.7s |
| reset-grid-keep-viewers | combine | pass | pass | pass | pass | **PASS** | 5.6s |
| sort-class-then-activity | combine | pass | pass | pass | pass | **PASS** | 4.0s |
| trap-sortbycolumns | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 8.8s |
| trap-gridcol-name-rename | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 7.2s |
| trap-reset-layout-ambiguous | wrong-tool-trap | pass | pass | pass | pass | **PASS** | 3.7s |
| reset-all-customization | reset-or-cleanup | pass | pass | pass | pass | **PASS** | 7.3s |
| clear-everything-fresh-start | reset-or-cleanup | pass | fail | fail | pass | **FAIL** | 6.9s |

## Failed prompts (deep)

### clear-everything-fresh-start — reset-or-cleanup

**Prompt:** Clear everything and start fresh
**Intent:** 'Everything' arguably includes added viewers. Pass if either a composed reset is performed, or view.resetLayout() is used (platform's blunt instrument for exactly this). Only fail if Claude does nothing or only changes one aspect.

**Rubric:**

- path: PASS
- helpers: FAIL (`view\.grid\.sort\(\s*\[\s*\]|setColumnsWidthType|setDisabled|\.visible\s*=\s*true|view\.resetLayout`)
- symbols: FAIL (`sort|reset|visible|clear`)
- forbidden: PASS

**Hypothesis:** Claude bypassed the expected helper(s) (`view\.grid\.sort\(\s*\[\s*\]|setColumnsWidthType|setDisabled|\.visible\s*=\s*true|view\.resetLayout`). It may have reached for the raw js-api equivalent — flag this anti-pattern more loudly in the skill or hide the raw helper. Expected symbol(s) absent (`sort|reset|visible|clear`). Either Claude used the wrong constant or chose a different valid path — review the response and decide whether to tighten the rubric or the skill.

**Emitted exec blocks:**

```javascript
grok.shell.closeAll();
```

<details><summary>Full response text</summary>

```
```datagrok-exec
grok.shell.closeAll();
```

Done! All views and tables have been closed. You're starting with a clean slate.
```

</details>

## Suggested skill improvements

- Expected helper `view\.grid\.sort\(\s*\[\s*\]|setColumnsWidthType|setDisabled|\.visible\s*=\s*true|view\.resetLayout` not reached for: clear-everything-fresh-start. Consider tightening the skill's routing table or example density for this case.

## Anomalies

- None.
