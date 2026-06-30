---
phase: 15
plan: 03
status: complete
type: implementation
---

# Plan 15-03 Summary — assertPublishedShape Round-Trip Gate

Built `src/publishing/assert-published-shape.ts`, the load-bearing round-trip gate consumed by BOTH the Plan 04 orchestrator (step 8 non-negotiable gate) and the Plan 08 regression suite. Single source of truth: if the contract changes, it changes here.

## Exports

- `FORMULA_LINE_ASSERTION_PREFIX = 'FORMULA_LINE_MISSING:'` — load-bearing prefix Plan 04's catch block matches on to invoke the W-7 self-healing recovery hook before rolling back
- `PublishedShapeContract` interface — 7 fields: expectedName, expectedProjectName, expectedAllowlist, expectedMeta, expectVolcano, expectEnrichment, expectFormulaLines
- `assertPublishedShape(project, contract): Promise<void>` — reopens via `dapi.projects.find(id).open()`, asserts 10 invariants (with 8a being the formula-line specialization), throws with a specific message on first violation

## Assertions

| # | Asserts | Throws on | Notes |
|---|---------|-----------|-------|
| 1 | `reDf.name === contract.expectedName` | got/expected mismatch | A5 — spike confirmed survival |
| 2 | `reopened.name === contract.expectedProjectName` | got/expected mismatch | PUB-09 |
| 3 | every `contract.expectedAllowlist[i]` is a column on `reDf` | column missing | one throw per missing col |
| 4 | visible (non-`_meta_*`) column count matches allowlist | unexpected extras / missing | enumerates both lists in the message |
| 5 | `findColumn(PROTEIN_ID).semType === SEMTYPE.PROTEIN_ID` | semType stripped | hint: re-assign in trim-dataframe Step E |
| 6 | `getPublishedMetadata(reDf)` returns non-null + matches contract on target / deMethod / fcThreshold (float eps 1e-9) / pThreshold / version / publishId / includesEnrichment / publishedBy | any field mismatch | column-FIRST tag-SECOND via Plan 01 helper |
| 7 | every key in `META_COLUMNS` is a column on `reDf` | missing belt-and-braces col | PUB-11 |
| 8 | volcano scatter plot present in `tv.viewers` | absent after 5000ms `awaitCheck` | hint: post-open `recomputeVolcano` |
| 8a | `df.meta.formulaLines.items` OR `volcano.getOptions().look.formulaLines` OR `volcano.props.formulaLines` contains a line referencing `fcThreshold` AND a line referencing `-log10(pThreshold)` (or `pThreshold` substring) | missing line(s) | throws with `FORMULA_LINE_ASSERTION_PREFIX` so Plan 04 step 8 detects and self-heals (B-2 / W-7) |
| 9 | `grok.shell.tables` has ≥ 2 DataFrames AND protein DF (`_meta_published_includes_enrichment === 'true'`) AND enrichment DF (`proteomics.enrichment` tag) | any missing | only runs when `contract.expectEnrichment` |
| 10 | `reopened.options[PUBLISHED_TAGS.SUPERSEDES]` and `[SUPERSEDED_BY]` match contract | got/expected mismatch | only runs when contract values non-null |

## How spike output drove the design

- Spike A6 confirmed scatter-plot survival with `xColumnName='log2FC'` and `yColumnName='adj.p-value'`. Assertion 8 verifies presence only (props vary across volcano variants); precise prop assertions are deferred to Plan 08's regression suite which controls the input shape.
- Spike showed all 4 set semTypes survived Project round-trip — but Assertion 5 still checks PROTEIN_ID explicitly because `trim-dataframe.ts` cloning is a separate boundary the spike did not exercise.
- Spike A1 confirmed `permissions.get` returns `["view", "edit"]` keys with empty body when no grants exist. Assertion 10 reads from `reopened.options[...]` not `permissions.get(...)` per Plan 04's grant-on-project decision in `15-00-SUMMARY` (Space-inheritance is NOT visible via `permissions.get(project)`).
- Spike Pitfall 3 result: all 14 tags survived. Assertion 6 reads via `getPublishedMetadata` which falls back column-first then tag, so even if a future serializer change strips tags, the assertion still passes via belt-and-braces columns.

## Single source of truth (Plan 04 + Plan 08)

Both Plan 04's `publishAnalysis` orchestrator step 8 and Plan 08's regression suite will:

```ts
import {assertPublishedShape, FORMULA_LINE_ASSERTION_PREFIX} from './assert-published-shape';
```

Plan 04 step 8 catch block matches on `err.message.startsWith(FORMULA_LINE_ASSERTION_PREFIX)` to trigger the post-open formula-line recovery (W-7).

## Verification

Project-wide `tsc --noEmit` passes clean. The file has 10 explicit `throw new Error(...)` sites covering every assertion in the contract. Each throw includes the failed-field name + observed + expected so failure diagnosis is single-step.

## Output

`src/publishing/assert-published-shape.ts` — 209 lines, type-checks under strict mode.
