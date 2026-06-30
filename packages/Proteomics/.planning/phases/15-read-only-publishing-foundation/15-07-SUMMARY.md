---
phase: 15
plan: 07
status: complete
type: implementation
---

# Plan 15-07 Summary — Package Wireup

Made the four minimum touch points to integrate the publishing primitives + reviewer panel + post-open recovery into the running package.

## Files changed

| File | Change |
|------|--------|
| `src/publishing/post-open-recovery.ts` | NEW — `recoverPublishedProject(df)` + reuses existing `wireEnrichmentToVolcano` from `src/viewers/enrichment-viewers.ts` |
| `src/package.ts` | +5 imports; +1 menu handler `shareAnalysisForReview`; +1 panel decorator `publishedAnalysisPanelWidget`; +1 autostart hook `recoverPublishedProjectsOnStartup` |
| `src/package-test.ts` | +1 import line for the Plan 08 regression test file |

## src/publishing/post-open-recovery.ts

Two-part self-healing for reopened published Projects:

1. **B-2 / W-7 — formula-line recovery.** Finds the TableView showing the DataFrame, locates the volcano scatter plot in it, checks whether formula lines are already present (via `df.meta.formulaLines.items` OR `volcano.getOptions().look.formulaLines`), and if absent reads FC and p thresholds from `proteomics.published_fc_threshold` / `proteomics.published_p_threshold` tags and re-applies via `applyVolcanoFormulaLines` (Plan 04 export).
2. **W-5 — enrichment subscription re-wire.** When a published enrichment DataFrame is in `grok.shell.tables` alongside the protein DataFrame (matching `proteomics.published === 'true'`), calls the existing `wireEnrichmentToVolcano` helper from `src/viewers/enrichment-viewers.ts` (no duplication — same export, single source of truth). Sentinel tag `proteomics.enrichment_wired === 'true'` prevents double subscription on repeat opens.

Both parts are wrapped in try/catch so per-DataFrame failures do not block recovery of other published DataFrames.

## src/package.ts edits

**Imports** added after `panels/uniprot-panel`:
```ts
import {publishedAnalysisPanel} from './panels/published-analysis-panel';
import {showShareForReviewDialog} from './publishing/share-dialog';
import {recoverPublishedProject} from './publishing/post-open-recovery';
import {isPublished} from './publishing/publish-state';
```

**Menu handler** `shareAnalysisForReview` registered as `Proteomics | Share | Share Analysis for Review...`:
- `requireDifferentialExpression` gate (CONTEXT.md Deferred Ideas item 11)
- `await showShareForReviewDialog(df)` per W-6 (the dialog is async because it pre-fetches `dapi.groups.list()`)

**Panel decorator** `publishedAnalysisPanelWidget` registered alongside `uniprotPanelWidget`:
- Name: `Proteomics | Shared Analysis` (Pitfall 14 — "Shared" not "Published")
- `semType=Proteomics-ProteinId` parameter filter
- Returns `publishedAnalysisPanel(proteinId)` — which adds its own first-line `isPublished` defense-in-depth

**Autostart hook** `recoverPublishedProjectsOnStartup` registered with `tags: ['autostart']` + `meta: {autostartImmediate: 'true'}` (string per the platform's meta-value typing):
- Tries `grok.events.onProjectOpened` first, falls back to `onCurrentProjectChanged`, finally `onViewAdded`
- On each event: 200ms settle, scan `grok.shell.tables` for `isPublished(df)`, call `recoverPublishedProject(df)` per published DF
- Per-DataFrame errors surface as `shell.warning` and continue to the next DF — never blocks startup

## src/package-test.ts edit

Added `import './tests/publish-roundtrip';` after the existing `import './tests/publish-spike';`. Plan 08 creates the referenced file; webpack will fail to bundle the test entry point until Plan 08 lands.

## SEMTYPE.DIRECTION decision

Per CONTEXT.md Claude's discretion + RESEARCH "NO new SEMTYPEs anticipated": no change to `src/utils/proteomics-types.ts` or `detectors.js`. The `direction` column stays string-typed and is matched by name hint in `findColumn(df, '', ['direction', 'regulation', 'up_down'])` (Plan 02's `trimForPublish`). v1.3 volcano consumes `direction` by name; no typed lookup needed.

## Verification

`tsc --noEmit -p tsconfig.json` passes clean for all source files. Webpack build fails on the deliberately-broken test import for `./tests/publish-roundtrip` — this is expected; Plan 08 lands the file and the build comes back clean.

`requireDifferentialExpression` is still the gate for behaviors that need DE-complete state.

## Output

- `src/publishing/post-open-recovery.ts` — 105 lines (NEW)
- `src/package.ts` — +66 lines (4 imports, 1 menu handler, 1 panel decorator, 1 autostart hook)
- `src/package-test.ts` — +1 import line
