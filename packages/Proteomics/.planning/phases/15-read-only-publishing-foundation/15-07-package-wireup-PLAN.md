---
phase: 15-read-only-publishing-foundation
plan: 07
type: execute
wave: 5
depends_on: ["15-04", "15-05", "15-06"]
files_modified:
  - src/package.ts
  - src/package-test.ts
  - src/publishing/post-open-recovery.ts
autonomous: true
requirements: [PUB-04, PUB-06, PUB-07, PUB-12]
must_haves:
  truths:
    - "Menu entry `Proteomics | Share | Share Analysis for Review...` appears in the top-menu and dispatches to `showShareForReviewDialog(df)` via `await` (W-6 — showShareForReviewDialog is async)"
    - "Menu handler gates on `requireDifferentialExpression(df)` — disabled/warning when active DF lacks `proteomics.de_complete === 'true'`"
    - "Panel `Proteomics | Shared Analysis` registered via `@grok.decorators.panel` with `semType=Proteomics-ProteinId` parameter filter"
    - "Post-open recovery hook (B-2 / W-7) registered as `@grok.decorators.func` with `meta.role: autostart` — fires on platform startup, subscribes to grok.events.onProjectOpened (or equivalent), detects `proteomics.published === 'true'` DataFrames, re-applies volcano formula lines from `proteomics.published_fc_threshold` / `proteomics.published_p_threshold` tags via the exported `applyVolcanoFormulaLines` helper"
    - "Same post-open recovery hook ALSO re-establishes the cross-DF enrichment-to-protein highlight subscription via `wireEnrichmentToVolcano` (W-5) when both DataFrames are present in `grok.shell.tables` with matching `proteomics.published` tags"
    - "`src/package-test.ts` imports `./tests/publish-roundtrip` so `grok test --category Publishing` resolves (Plan 08 creates the test file)"
  artifacts:
    - path: "src/publishing/post-open-recovery.ts"
      provides: "Reopened-Project recovery hook: re-applies formula lines from tags + re-wires enrichment cross-DF subscription"
      exports: ["recoverPublishedProject", "wireEnrichmentToVolcano"]
    - path: "src/package.ts"
      provides: "Share menu handler + published-analysis panel decorator + autostart registration of recovery hook"
      contains: "shareAnalysisForReview, publishedAnalysisPanelWidget, recoverPublishedProjectOnOpen"
    - path: "src/package-test.ts"
      provides: "Publish roundtrip test registration"
      contains: "import './tests/publish-roundtrip';"
  key_links:
    - from: "src/package.ts"
      to: "src/publishing/share-dialog.ts"
      via: "showShareForReviewDialog import"
      pattern: "import.*share-dialog"
    - from: "src/package.ts"
      to: "src/panels/published-analysis-panel.ts"
      via: "publishedAnalysisPanel import"
      pattern: "import.*published-analysis-panel"
    - from: "src/package.ts"
      to: "src/publishing/post-open-recovery.ts"
      via: "recoverPublishedProject import + autostart registration"
      pattern: "recoverPublishedProject"
---

<objective>
Wire the Wave 3 surfaces (share dialog + reviewer panel) into the package entry point so they're reachable from the running platform. Single menu item + single panel decorator + single test import line. This is the integration moment where the publishing primitives become available to users.

**Revision B-2 / W-5 / W-7 additions:**
- Create new file `src/publishing/post-open-recovery.ts` exporting:
  - `recoverPublishedProject(df: DG.DataFrame): Promise<void>` — given a reopened DataFrame that has `proteomics.published === 'true'`, re-apply volcano formula lines from the threshold tags AND re-establish the cross-DF enrichment subscription when present.
  - `wireEnrichmentToVolcano(enrichDf: DG.DataFrame, proteinDf: DG.DataFrame): void` — the cross-DF current-row subscription pattern from v1.2 `src/viewers/enrichment-viewers.ts`, lifted into a reusable helper.
- Register the recovery hook as an `@grok.decorators.func` with `meta.role: 'autostart'` (or `'init'`) in `src/package.ts` so it fires on platform startup and subscribes to project-open events. On every project open: scan `grok.shell.tables` for published DataFrames and call `recoverPublishedProject(df)` on each.

Purpose: Plans 01-06 built the publishing primitives + reviewer-side panel in isolation; they don't appear in the UI yet, and the reopen-recovery path is not wired. This plan adds the four minimum touch points: (1) the `@grok.decorators.func` for the Share menu, (2) the `@grok.decorators.panel` for the reviewer audit context, (3) the autostart-registered post-open recovery hook (B-2 / W-5), (4) the test-import line so Plan 08's regression test gets picked up by `grok test`.

Output: new file `src/publishing/post-open-recovery.ts` + edits to `src/package.ts` (3 imports + menu handler + panel decorator + autostart hook) + 1 edit to `src/package-test.ts` (import line).
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/package.ts
@packages/Proteomics/src/package-test.ts
@packages/Proteomics/src/publishing/share-dialog.ts
@packages/Proteomics/src/publishing/publish-project.ts
@packages/Proteomics/src/publishing/publish-state.ts
@packages/Proteomics/src/panels/published-analysis-panel.ts
@packages/Proteomics/src/viewers/volcano.ts
@packages/Proteomics/src/viewers/enrichment-viewers.ts
@packages/Proteomics/src/analysis/differential-expression.ts

<interfaces>
<!-- From src/package.ts existing surface -->

From src/analysis/differential-expression.ts (existing):
  requireDifferentialExpression(df: DG.DataFrame, msg: string): boolean
  // returns true if df.getTag('proteomics.de_complete') === 'true'; else shows warning(msg) and returns false

From src/package.ts decorator usage:
  @grok.decorators.func({'top-menu': 'Proteomics | Group | Item Label...'})
  static async handler(): Promise<void> { ... }

  @grok.decorators.panel({name: 'Proteomics | Foo', description: '...', meta: {role: 'widgets'}})
  static panelWidget(@grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string): DG.Widget { ... }

  @grok.decorators.func({tags: ['autostart'], meta: {autostartImmediate: true}})
  static async autostartHandler(): Promise<void> { ... }

From Plan 05:
  showShareForReviewDialog(df: DG.DataFrame): Promise<void>   // W-6: async

From Plan 06:
  publishedAnalysisPanel(proteinId: string): DG.Widget

From Plan 04:
  applyVolcanoFormulaLines(viewer: DG.Viewer, fcThreshold: number, pThreshold: number): void

From Plan 01:
  isPublished(df: DG.DataFrame): boolean
  PUBLISHED_TAGS, DE_COMPLETE_TAG

From src/viewers/enrichment-viewers.ts (existing — v1.2 precedent for cross-DF subscription):
  (look for the onCurrentRowChanged subscription that highlights matching proteins in the protein DF
   from a selected enrichment term row — extract into the new wireEnrichmentToVolcano helper)

From datagrok-api:
  grok.events.onProjectOpened: rxjs.Observable<DG.Project>    // or equivalent — verify against platform
  grok.shell.tables: DG.DataFrame[]
  grok.shell.tv: DG.TableView | null
  v.type: string; DG.VIEWER.SCATTER_PLOT = 'Scatter plot'
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Create post-open-recovery.ts with recoverPublishedProject + wireEnrichmentToVolcano</name>
  <files>src/publishing/post-open-recovery.ts</files>
  <read_first>
    - @packages/Proteomics/src/publishing/publish-state.ts (PUBLISHED_TAGS — PUBLISHED, PUBLISHED_FC_THRESHOLD, PUBLISHED_P_THRESHOLD)
    - @packages/Proteomics/src/publishing/publish-project.ts (applyVolcanoFormulaLines export)
    - @packages/Proteomics/src/viewers/enrichment-viewers.ts (full file — v1.2 cross-DF onCurrentRowChanged precedent; this is the pattern to extract into wireEnrichmentToVolcano)
    - @packages/Proteomics/src/viewers/volcano.ts (volcano viewer type identification)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-05 enrichment carry; PUB-06 SC-2 formula lines mandate)
    - @packages/Proteomics/CLAUDE.md (rxjs subscription cleanup convention)
  </read_first>
  <action>
Create `src/publishing/post-open-recovery.ts`. Standard imports plus `isPublished`, `PUBLISHED_TAGS` from `./publish-state`, `applyVolcanoFormulaLines` from `./publish-project`.

**Export `recoverPublishedProject(df: DG.DataFrame): Promise<void>`** — the load-bearing self-healing function. Body:

```
export async function recoverPublishedProject(df: DG.DataFrame): Promise<void> {
  // Guard: only act on published DFs.
  if (!isPublished(df)) return;
  
  // PART 1 — Re-apply formula lines on the volcano (B-2 / W-7 belt-and-braces #3)
  // Find the TableView showing this DF, find the volcano viewer in it.
  const tv = grok.shell.tableViews.find((view) => view.dataFrame === df);
  if (tv) {
    const volcano = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
    if (volcano) {
      // Read thresholds from belt-and-braces tags (which Plan 02 setPublishedTags wrote, and which survive
      // serializer even when look config doesn't per Phase 13 e527d07ba1 evidence).
      const fcRaw = df.getTag(PUBLISHED_TAGS.PUBLISHED_FC_THRESHOLD);
      const pRaw = df.getTag(PUBLISHED_TAGS.PUBLISHED_P_THRESHOLD);
      if (fcRaw && pRaw) {
        const fc = parseFloat(fcRaw);
        const p = parseFloat(pRaw);
        if (!isNaN(fc) && !isNaN(p)) {
          // Check whether formula lines are already present — if so, no-op (idempotent).
          const opts = volcano.getOptions?.() as any;
          const existingLines = opts?.look?.formulaLines;
          if (!existingLines || existingLines === '[]' || existingLines === '') {
            applyVolcanoFormulaLines(volcano, fc, p);
          }
        }
      }
    }
  }
  
  // PART 2 — Re-establish enrichment cross-DF subscription (W-5)
  // When the published Project carries BOTH protein DF AND enrichment DF, re-wire the highlight subscription.
  const enrichDf = grok.shell.tables.find((t) => 
    t !== df && 
    t.getTag('proteomics.enrichment') === 'true' && 
    t.getTag(PUBLISHED_TAGS.PUBLISHED) === 'true'
  );
  if (enrichDf) {
    wireEnrichmentToVolcano(enrichDf, df);
  }
}
```

**Export `wireEnrichmentToVolcano(enrichDf: DG.DataFrame, proteinDf: DG.DataFrame): void`** — extracted from the v1.2 `enrichment-viewers.ts` precedent. Reads the existing subscription pattern (subscribe to `enrichDf.onCurrentRowChanged`, parse the intersection-genes string from the current row, set selection on `proteinDf` matching those gene symbols / protein IDs).

To avoid duplicate subscriptions on re-trigger (recovery hook may fire multiple times if the user reopens the project repeatedly in one session), store a sentinel tag `proteomics.enrichment_wired === 'true'` on the enrichment DF after wiring; the function early-returns if the sentinel is already present.

```
export function wireEnrichmentToVolcano(enrichDf: DG.DataFrame, proteinDf: DG.DataFrame): void {
  // Sentinel guard against duplicate subscriptions
  if (enrichDf.getTag('proteomics.enrichment_wired') === 'true') return;
  
  // Lift the existing subscription pattern from src/viewers/enrichment-viewers.ts
  enrichDf.onCurrentRowChanged.subscribe(() => {
    const currentIdx = enrichDf.currentRowIdx;
    if (currentIdx < 0) return;
    
    // Extract intersection-genes from the current enrichment row
    // (column name: 'Intersection' per the enrichment DF column convention)
    const intersectionCol = enrichDf.col('Intersection');
    if (!intersectionCol) return;
    const genesStr = intersectionCol.get(currentIdx) as string;
    if (!genesStr) return;
    
    // Parse comma-separated gene list and highlight matching proteins in the protein DF
    const genes = genesStr.split(',').map((g) => g.trim()).filter(Boolean);
    const geneCol = proteinDf.col('Gene Name') || proteinDf.col('Gene names');
    if (!geneCol) return;
    
    proteinDf.selection.init((rowIdx) => {
      const g = geneCol.get(rowIdx) as string;
      return !!g && genes.includes(g);
    });
  });
  
  enrichDf.setTag('proteomics.enrichment_wired', 'true');
}
```

**Note:** the exact column name / pattern depends on what `src/viewers/enrichment-viewers.ts` does today. Executor should READ that file first and lift the actual pattern — the snippet above is a structural sketch.

JSDoc: "Post-open recovery hook. Self-heals reopened published Projects by re-applying volcano formula lines from belt-and-braces tags (B-2 / W-7) and re-establishing the enrichment cross-DF highlight subscription (W-5). Idempotent — safe to call multiple times on the same DataFrame."
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/post-open-recovery.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export async function recoverPublishedProject" src/publishing/post-open-recovery.ts | grep -q '^1$' &amp;&amp; grep -c "export function wireEnrichmentToVolcano" src/publishing/post-open-recovery.ts | grep -q '^1$' &amp;&amp; grep -c "applyVolcanoFormulaLines" src/publishing/post-open-recovery.ts | grep -qv '^0$' &amp;&amp; grep -c "PUBLISHED_FC_THRESHOLD" src/publishing/post-open-recovery.ts | grep -qv '^0$' &amp;&amp; grep -c "onCurrentRowChanged" src/publishing/post-open-recovery.ts | grep -qv '^0$' &amp;&amp; grep -c "enrichment_wired" src/publishing/post-open-recovery.ts | grep -qv '^0$'</automated>
  </verify>
  <done>New file `src/publishing/post-open-recovery.ts` exports `recoverPublishedProject` + `wireEnrichmentToVolcano`; uses `isPublished` guard; re-applies formula lines via `applyVolcanoFormulaLines` from tags (idempotent — no-ops if already present); re-establishes cross-DF subscription with sentinel-tag duplicate-guard; strict TypeScript compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 2: Add Share menu handler (await async) + panel decorator + autostart recovery hook in src/package.ts</name>
  <files>src/package.ts</files>
  <read_first>
    - @packages/Proteomics/src/package.ts (lines 359-371 — DE menu handler precedent; lines 497-512 — heatmap precondition-gate precedent; search for existing autostart functions to find the right slot)
    - @packages/Proteomics/src/analysis/differential-expression.ts (`requireDifferentialExpression` signature + behavior)
    - @packages/Proteomics/src/publishing/publish-state.ts (PUBLISHED_TAGS, isPublished)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 9 — menu wireup analog + notes)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (Deferred Ideas item 11 — publishing non-DE-complete DFs is out of scope; menu gates on de_complete)
    - @packages/Proteomics/src/publishing/post-open-recovery.ts (just-built recoverPublishedProject signature)
  </read_first>
  <action>
Open `src/package.ts`. Make 4 minimal edits:

**Edit 1 — add imports:**
- `import {showShareForReviewDialog} from './publishing/share-dialog';`
- `import {publishedAnalysisPanel} from './panels/published-analysis-panel';`
- `import {recoverPublishedProject} from './publishing/post-open-recovery';`
- `import {isPublished} from './publishing/publish-state';`

Place these in the existing import block; group them after the `panels/uniprot-panel` import to keep similar surfaces adjacent.

**Edit 2 — add Share menu handler (model after the DE menu handler at lines 359-371; W-6: await async dialog):**

Inside the `PackageFunctions` class, add a new static async method:

```
@grok.decorators.func({'top-menu': 'Proteomics | Share | Share Analysis for Review...'})
static async shareAnalysisForReview(): Promise<void> {
  const tv = grok.shell.tv;
  const df = tv?.dataFrame;
  if (!tv || !df) { grok.shell.warning('No table open'); return; }
  if (!requireDifferentialExpression(df,
    'Run Differential Expression first (Proteomics | Analyze | Differential Expression)')) return;
  await showShareForReviewDialog(df);  // W-6: showShareForReviewDialog is async
}
```

Location: insert AFTER the heatmap handler (`showHeatmap`) and BEFORE the panel-decorator block — keep all menu handlers contiguous. The new `Share` submenu has no other entries (CONTEXT.md D-04: "No separate 'Update Share' menu item — keeps menu surface flat").

Per CLAUDE.md naming convention: handler method is `shareAnalysisForReview` (camelCase); menu path ends in `...` per dialog-suffix convention.

**Edit 3 — add Published Analysis panel decorator (model after the UniProt panel at lines 630-639):**

```
@grok.decorators.panel({
  name: 'Proteomics | Shared Analysis',
  description: 'Audit context for a shared analysis snapshot',
  meta: {role: 'widgets'},
})
static publishedAnalysisPanelWidget(
  @grok.decorators.param({options: {semType: 'Proteomics-ProteinId'}}) proteinId: string,
): DG.Widget {
  return publishedAnalysisPanel(proteinId);
}
```

Location: insert alongside the existing `uniprotPanelWidget` panel decorator (lines 630-639). Per Pitfall 14: panel name is "Shared Analysis" not "Published Analysis" (biologist-jargon-natural).

**Edit 4 — register autostart post-open recovery hook (B-2 / W-5 / W-7):**

```
@grok.decorators.func({tags: ['autostart'], meta: {autostartImmediate: true}})
static async recoverPublishedProjectsOnStartup(): Promise<void> {
  // Subscribe to project-open events; on each open, scan shell tables for published DFs and self-heal.
  // The exact event API depends on Datagrok release — try grok.events.onProjectOpened first;
  // if not available, fall back to grok.events.onCurrentProjectChanged or onViewAdded with isPublished filter.
  const tryEvent = (grok.events as any).onProjectOpened ?? (grok.events as any).onCurrentProjectChanged ?? null;
  if (tryEvent && typeof tryEvent.subscribe === 'function') {
    tryEvent.subscribe(async () => {
      // Allow shell.tables to populate after the project open event fires
      await new Promise((r) => setTimeout(r, 200));
      for (const df of grok.shell.tables) {
        if (isPublished(df)) {
          try { await recoverPublishedProject(df); } catch (e: any) {
            grok.shell.warning('Could not auto-recover shared analysis on open: ' + (e?.message ?? e));
          }
        }
      }
    });
  } else {
    // Defensive fallback: subscribe to grok.events.onViewAdded (always available) and check the new view's DF.
    grok.events.onViewAdded.subscribe(async (view: DG.View) => {
      if (view instanceof DG.TableView) {
        const df = view.dataFrame;
        if (df && isPublished(df)) {
          try { await recoverPublishedProject(df); } catch (e: any) {
            grok.shell.warning('Could not auto-recover shared analysis on open: ' + (e?.message ?? e));
          }
        }
      }
    });
  }
}
```

Location: insert near the bottom of the class, with other autostart/init functions if any exist. Per CLAUDE.md function-role table: `#autostart` tag registers it as a startup function.

After all four edits, verify:
- No duplicate `@grok.decorators.func` block for `Proteomics | Share | ...`
- No duplicate `@grok.decorators.panel` block for `Proteomics | Shared Analysis`
- No duplicate `recoverPublishedProjectsOnStartup` autostart function
- All imports resolve (Plan 01-06 + the new post-open-recovery.ts files exist)

Per CLAUDE.md: never edit `.g.ts` or `.api.g.ts` files manually. After this plan completes, the build pipeline runs `grok api && grok check` which regenerates `src/package.g.ts` — that regen step is part of `npm run build` (NOT this plan's responsibility, but the executor should run `npm run build` after edits to confirm registration).
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; grep -c "import {showShareForReviewDialog}" src/package.ts | grep -q '^1$' &amp;&amp; grep -c "import {publishedAnalysisPanel}" src/package.ts | grep -q '^1$' &amp;&amp; grep -c "import {recoverPublishedProject}" src/package.ts | grep -q '^1$' &amp;&amp; grep -c "shareAnalysisForReview" src/package.ts | grep -qv '^0$' &amp;&amp; grep -c "await showShareForReviewDialog" src/package.ts | grep -qv '^0$' &amp;&amp; grep -c "Share Analysis for Review" src/package.ts | grep -qv '^0$' &amp;&amp; grep -c "publishedAnalysisPanelWidget" src/package.ts | grep -qv '^0$' &amp;&amp; grep -c "recoverPublishedProjectsOnStartup\\|recoverPublishedProject" src/package.ts | grep -qv '^0$' &amp;&amp; grep -c "requireDifferentialExpression" src/package.ts | grep -qv '^0$' &amp;&amp; npx tsc --noEmit src/package.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; }</automated>
  </verify>
  <done>`src/package.ts` contains imports for `showShareForReviewDialog`, `publishedAnalysisPanel`, `recoverPublishedProject`, and `isPublished`; menu handler `shareAnalysisForReview` registered with `Proteomics | Share | Share Analysis for Review...` path; uses `await showShareForReviewDialog(df)` (W-6); `requireDifferentialExpression` gate present; panel decorator `publishedAnalysisPanelWidget` registered with `semType=Proteomics-ProteinId`; autostart hook `recoverPublishedProjectsOnStartup` registered that fires on project-open events and self-heals via `recoverPublishedProject` (B-2 / W-5 / W-7); strict TypeScript compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 3: Register publish-roundtrip test in src/package-test.ts (conditional SEMTYPE.DIRECTION decision)</name>
  <files>src/package-test.ts</files>
  <read_first>
    - @packages/Proteomics/src/package-test.ts (existing import block lines 4-20 + Plan 00's edit adding publish-spike import)
    - @packages/Proteomics/src/utils/proteomics-types.ts (current SEMTYPE list — check whether DIRECTION already exists)
    - @packages/Proteomics/.claude/feedback_grok_test_skipbuild_stale.md (memory: rebuild after adding new test files)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (Claude's discretion: direction column representation — string by default, no SEMTYPE needed unless v1.3 volcano consumes a typed version)
  </read_first>
  <action>
Open `src/package-test.ts`. Append `import './tests/publish-roundtrip';` at the end of the import block (after the spike import added by Plan 00). Match the existing single-quote + semicolon convention.

**Conditional decision — SEMTYPE.DIRECTION:**

Per CONTEXT.md Claude's discretion + RESEARCH §"Recommended Project Structure" line 164: default is NO new SEMTYPEs. The direction column stays string-typed.

Open `src/utils/proteomics-types.ts` to check:
- If `SEMTYPE.DIRECTION` already exists: no change needed.
- If `SEMTYPE.DIRECTION` does NOT exist: check whether `src/viewers/volcano.ts` references `SEMTYPE.DIRECTION`. If not, leave both files unchanged (default per RESEARCH).

If a determination is made to ADD `SEMTYPE.DIRECTION`:
- Update `src/utils/proteomics-types.ts` adding `DIRECTION: 'Proteomics-Direction',` to the `SEMTYPE` object.
- Update `detectors.js` (root-level, NOT under `src/`) adding a `detectDirection(col)` block per PATTERNS.md Section 11 analog:
  ```
  //meta.role: semTypeDetector
  //input: column col
  //output: string semType
  detectDirection(col) {
    if (col.type !== DG.TYPE.STRING) return null;
    const name = col.name.toLowerCase();
    if (name === 'direction' || name === 'regulation') {
      col.semType = 'Proteomics-Direction';
      return col.semType;
    }
    return null;
  }
  ```
- Update `src/publishing/trim-dataframe.ts` Step E (Plan 02) to re-assign `SEMTYPE.DIRECTION` on the cloned direction column.

DEFAULT (per RESEARCH "NO new SEMTYPEs anticipated"): leave both files unchanged. Document the decision in this plan's SUMMARY.

After edits, per memory `feedback_grok_test_skipbuild_stale.md`: rebuild the test bundle (`npm run build`, NOT `grok test --skip-build`) before invoking the Publishing tests for the first time.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; grep -c "import './tests/publish-roundtrip';" src/package-test.ts | grep -q '^1$'</automated>
  </verify>
  <done>`src/package-test.ts` contains exactly one `import './tests/publish-roundtrip';` line; SEMTYPE.DIRECTION decision documented (default: no change to `proteomics-types.ts` or `detectors.js`).</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| menu click -> showShareForReviewDialog | Menu handler gates on `requireDifferentialExpression`; dispatches to Plan 05 dialog (now async — handler awaits) |
| panel decorator -> publishedAnalysisPanel(proteinId) | Platform invokes panel with semType-filtered protein ID; panel's first-line `isPublished` guard provides defense-in-depth |
| autostart hook -> reopened published DataFrames | Hook subscribes to project-open events; runs self-heal only on DFs that pass `isPublished` check; idempotent re-application |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-SP-10 | Elevation of privilege (publishing without DE complete) | shareAnalysisForReview menu handler | mitigate | `requireDifferentialExpression(df, msg)` gate prevents publish from running on annotated-only / normalized-only DataFrames; CONTEXT.md "Deferred Ideas" item 11 confirms this is out of scope |
| T-15-SP-11 | Information disclosure (panel renders on non-published DF — defense-in-depth) | publishedAnalysisPanelWidget | mitigate | Panel decorator filters to `semType=Proteomics-ProteinId`; Plan 06 panel function adds first-line `isPublished` guard; double-layer defense |
| T-15-SP-12 | Information disclosure (formula lines stripped — recovery hook self-heals) | recoverPublishedProjectsOnStartup autostart hook | mitigate | Hook fires on every project open, scans shell tables for `isPublished(df)`, re-applies formula lines from belt-and-braces tags (B-2 belt-and-braces #3); idempotent — checks `existingLines` before re-apply |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- After `npm run build`: `grok api` regenerates `src/package.g.ts` with the new function/panel/autostart registrations
- `grok s functions list --filter "Proteomics:"` shows `shareAnalysisForReview`, `publishedAnalysisPanelWidget`, and `recoverPublishedProjectsOnStartup`
- Manual menu probe: open package, verify `Proteomics → Share → Share Analysis for Review...` appears
- Manual recovery probe: open a published Project; volcano formula lines should appear (either survived the round trip OR were re-applied by the hook within ~200ms of view materialization)
- `grok test --category Publishing-Spike` still works (Plan 00 spike import unaffected)
</verification>

<success_criteria>
- 4 import lines added to `src/package.ts`
- 1 menu handler `shareAnalysisForReview` registered, uses `await showShareForReviewDialog(df)` (W-6)
- 1 panel decorator `publishedAnalysisPanelWidget` registered
- 1 autostart hook `recoverPublishedProjectsOnStartup` registered (B-2 / W-5 / W-7)
- 1 new file `src/publishing/post-open-recovery.ts` with `recoverPublishedProject` + `wireEnrichmentToVolcano`
- 1 test import line added to `src/package-test.ts`
- SEMTYPE.DIRECTION decision documented (default: no change)
- `npm run build` succeeds without TypeScript errors
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-07-SUMMARY.md` when done with:
- Confirmation of the 4 minimal edits in src/package.ts + new file post-open-recovery.ts + 1 edit in src/package-test.ts
- SEMTYPE.DIRECTION decision (default: no change) recorded
- Build output confirmed clean
- Confirmation that the autostart recovery hook is wired via grok.events (project-open OR onViewAdded fallback) and idempotently re-applies formula lines + re-wires enrichment subscription
- Note for Plan 08 executor: `grok test --category Publishing` should now resolve once the roundtrip test file is created; Test 2 will exercise the wireEnrichmentToVolcano path
</output>
