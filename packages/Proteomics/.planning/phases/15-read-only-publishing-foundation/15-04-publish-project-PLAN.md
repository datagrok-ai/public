---
phase: 15-read-only-publishing-foundation
plan: 04
type: execute
wave: 3
depends_on: ["15-01", "15-02", "15-03"]
files_modified:
  - src/publishing/publish-project.ts
autonomous: true
requirements: [PUB-01, PUB-05, PUB-06, PUB-09, PUB-10, PUB-12]
must_haves:
  truths:
    - "Calling `publishAnalysis(df, opts)` produces a saved `DG.Project` named `Proteomics-Review-<slug>-v<N>-<YYYY-MM-DD>` containing the trimmed clone + optional enrichment DF"
    - "Reviewer group has VIEW-ONLY access on the per-target child Space; verify-and-rollback gate aborts publish (deletes Project, throws exact error string) if Edit slipped in via Space inheritance"
    - "`assertPublishedShape` is called inside `publishAnalysis` before returning success; failure rolls back via `dapi.projects.delete(project)`"
    - "Volcano formula lines for FC + p thresholds are re-applied AFTER the trim clone + TableView open BUT BEFORE the save (PUB-06 / SC-2 belt-and-braces #1), AND the threshold values are persisted in proteomics.published_fc_threshold / proteomics.published_p_threshold tags (PUB-06 / SC-2 belt-and-braces #2 so the post-open recovery hook can re-apply them if the serializer strips look/filter config)"
    - "Step 8 catch detects the FORMULA_LINE_ASSERTION_PREFIX failure mode SPECIFICALLY and invokes the post-open recovery hook once (re-applies formula lines from tags) before retrying the assertion — only if STILL failing does it roll back (W-7 self-healing path)"
    - "Republish creates a NEW Project with bidirectional supersede pointers (`proteomics.supersedes` on the new, `proteomics.superseded_by` on the prior); prior version is NOT deleted; supersede pointer is written to BOTH the prior Project's options AND the prior Project's DataFrame tag (W-8 belt-and-braces)"
    - "Source DataFrame is never mutated during the publish (Pitfall 1 guarantee)"
    - "Optional `opts.umbrellaName` override is supported (defaults to 'Proteomics-Reviews') so Plan 08 Test 5 can inject a throwaway umbrella to exercise the verify-and-rollback negative path"
  artifacts:
    - path: "src/publishing/publish-project.ts"
      provides: "publishAnalysis orchestrator + ensureUmbrellaSpace + ensurePerTargetChildSpace + verifyAndRollback + supersedeChain + applyVolcanoFormulaLines"
      exports: ["publishAnalysis", "applyVolcanoFormulaLines"]
  key_links:
    - from: "src/publishing/publish-project.ts"
      to: "src/publishing/trim-dataframe.ts"
      via: "trimForPublish + trimEnrichmentForPublish imports"
      pattern: "import.*trim-dataframe"
    - from: "src/publishing/publish-project.ts"
      to: "src/publishing/assert-published-shape.ts"
      via: "assertPublishedShape + FORMULA_LINE_ASSERTION_PREFIX imports + call before return"
      pattern: "assertPublishedShape\\("
    - from: "src/publishing/publish-project.ts"
      to: "grok.dapi.permissions"
      via: "grant + get + verify-and-rollback"
      pattern: "permissions\\.(grant|get)"
    - from: "src/publishing/publish-project.ts"
      to: "src/viewers/volcano.ts"
      via: "createVolcanoPlot import for B-2 re-render path"
      pattern: "createVolcanoPlot"
---

<objective>
Implement the orchestrator `publishAnalysis(df, opts)` that sequences the publish flow: trim → trim enrichment (opportunistic) → ensure umbrella Space → ensure per-target child Space → trim-side volcano re-creation + formula-line application (B-2) → save Project → grant View → verify-and-rollback → assertPublishedShape (with W-7 self-healing on FORMULA_LINE_MISSING) → write supersede chain (W-8 dual-write). Two non-negotiable gates (steps 7 and 8) each roll back via `dapi.projects.delete(project)` on failure with the exact error strings from CONTEXT.md D-03.

Purpose: this is the integration point where all the Wave 1 primitives compose into a single function the dialog handler (Plan 05) calls. Sequence is non-negotiable per RESEARCH §"System Architecture Diagram". Steps 7 (verify-and-rollback) and 8 (assert-published-shape) are the TWO non-negotiable gates per CONTEXT.md — neither can be skipped, both roll back on failure.

**Revision B-2 / W-7 addition:** PUB-06 / SC-2 requires formula lines to render on reviewer first paint. Phase 13 e527d07ba1 proves the serializer can strip look / filter config. Belt-and-braces mitigation:
  1. Persist `meta.fcThreshold` and `meta.pThreshold` as `proteomics.published_fc_threshold` / `proteomics.published_p_threshold` tags via `setPublishedTags` (already in Plan 01 PUBLISHED_TAGS) so a post-open recovery hook can re-apply them.
  2. Before the save (step 6.5), re-create the volcano on the trimmed DataFrame via the existing `createVolcanoPlot(frozen)` factory and apply formula lines for the two thresholds via `applyVolcanoFormulaLines(viewer, fcThreshold, pThreshold)`.
  3. After save + reopen (step 8), if `assertPublishedShape` throws with `FORMULA_LINE_ASSERTION_PREFIX`, the orchestrator invokes the post-open recovery hook ONCE (re-apply formula lines from the published tags, re-save the view) and retries the assertion. Only if STILL failing does it roll back.
The post-open recovery hook itself (an autostart-registered `@grok.decorators.func` that fires on Project open and self-heals reopened published Projects) is registered in Plan 07. This orchestrator OWNS the publish-side write of the threshold tags + the volcano re-render before save; the post-open hook OWNS the read-from-tag + re-apply path on reopen.

Output: `src/publishing/publish-project.ts` exporting `publishAnalysis(df, opts) -> Promise<DG.Project>` + `applyVolcanoFormulaLines(viewer, fcThreshold, pThreshold) -> void` (exported so Plan 07's post-open recovery hook can re-use it).
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md
@.planning/phases/15-read-only-publishing-foundation/15-01-SUMMARY.md
@.planning/phases/15-read-only-publishing-foundation/15-02-SUMMARY.md
@.planning/phases/15-read-only-publishing-foundation/15-03-SUMMARY.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/publishing/publish-state.ts
@packages/Proteomics/src/publishing/trim-dataframe.ts
@packages/Proteomics/src/publishing/assert-published-shape.ts
@packages/Proteomics/src/viewers/volcano.ts
@packages/ApiSamples/scripts/dapi/projects.js
@packages/ApiSamples/scripts/dapi/spaces.js
@packages/Bio/src/tests/projects-tests.ts

<interfaces>
<!-- Pulled from publishing primitives + Datagrok dapi -->

From Plan 01 (publish-state.ts):
  PublishOptions: { target, reviewerGroup, note, priorVersion, umbrellaName? }
  PublishedMetadata: { target, publishedAt, publishedBy, publishedByEmail, deMethod, fcThreshold, pThreshold, version, publishId, includesEnrichment, supersedes, supersededBy, ... }
  PUBLISHED_TAGS (includes PUBLISHED_FC_THRESHOLD, PUBLISHED_P_THRESHOLD, SUPERSEDED_BY, SUPERSEDES), META_COLUMNS
  slugifyTarget(raw): string
  findPriorShare(target, group): Promise<DG.Project | null>
  isPublished(df): boolean
  setPublishedTags(df, meta): void

From Plan 02 (trim-dataframe.ts):
  trimForPublish(source: DG.DataFrame, meta: PublishedMetadata): DG.DataFrame
  trimEnrichmentForPublish(enrichSource: DG.DataFrame, meta: PublishedMetadata): DG.DataFrame

From Plan 03 (assert-published-shape.ts):
  assertPublishedShape(project: DG.Project, contract: PublishedShapeContract): Promise<void>
  PublishedShapeContract: { expectedName, expectedProjectName, expectedAllowlist, expectedMeta, expectVolcano, expectEnrichment, expectFormulaLines }
  FORMULA_LINE_ASSERTION_PREFIX: 'FORMULA_LINE_MISSING:'

From src/viewers/volcano.ts:
  createVolcanoPlot(df: DG.DataFrame): DG.Viewer

From datagrok-api:
  DG.Project.create(): DG.Project
  project.name: string; project.addChild(entity); project.options: MapProxy
  grok.dapi.tables.uploadDataFrame(df): Promise<void>
  grok.dapi.tables.save(tableInfo): Promise<DG.TableInfo>
  grok.dapi.views.save(viewInfo): Promise<DG.ViewInfo>
  grok.dapi.projects.save(project): Promise<DG.Project>
  grok.dapi.projects.delete(project): Promise<void>
  grok.dapi.projects.find(id): Promise<DG.Project>
  grok.dapi.projects.filter(smart): HttpDataSource<DG.Project>
  grok.dapi.spaces.rootSpaceExists(name): Promise<boolean>
  grok.dapi.spaces.createRootSpace(name): Promise<DG.Project>      // Spaces ARE Projects with isSpace=true
  grok.dapi.spaces.id(id): SpaceClient
  spaceClient.subspaceExists(name): Promise<boolean>
  spaceClient.addSubspace(name): Promise<DG.Project>
  spaceClient.addEntity(entityId: string): Promise<void>
  spaceClient.children: SpaceChildrenClient { filter, list }
  grok.dapi.permissions.grant(entity, group, edit: boolean): Promise<void>
  grok.dapi.permissions.get(entity): Promise<{view: DG.Group[]; edit: DG.Group[]}>  // shape locked by spike A1
  grok.shell.user.friendlyName: string
  grok.shell.user.email: string | null
  grok.shell.tables: DG.DataFrame[]
  grok.shell.tv: DG.TableView | null     // active view (for layoutInfo)
  DG.TaskBarProgressIndicator.create(label): { update(msg, percent?), close() }
  viewer.setOptions(opts): void
  viewer.getOptions(): {look?: any; ...}
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Implement publishAnalysis steps 1-5 (trim + Spaces + save Project) + applyVolcanoFormulaLines helper</name>
  <files>src/publishing/publish-project.ts</files>
  <read_first>
    - @packages/ApiSamples/scripts/dapi/projects.js (canonical 3-line save shape)
    - @packages/Bio/src/tests/projects-tests.ts (lines 26-47 — multi-child save + reopen pattern)
    - @packages/ApiSamples/scripts/dapi/spaces.js (lines 6-67 — rootSpaceExists/createRootSpace + addSubspace + addEntity)
    - @packages/Proteomics/src/publishing/publish-state.ts (PublishOptions, PUBLISHED_TAGS including PUBLISHED_FC_THRESHOLD/PUBLISHED_P_THRESHOLD, slugifyTarget)
    - @packages/Proteomics/src/publishing/trim-dataframe.ts (trimForPublish, trimEnrichmentForPublish signatures)
    - @packages/Proteomics/src/viewers/volcano.ts (full file — createVolcanoPlot signature + how formula lines are added today)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 3 — orchestrator analogs A + B + C)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"System Architecture Diagram" — full sequence; §"Pattern 1" + §"Pattern 2")
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-03 nested Spaces — primary shape verified; PUB-06 SC-2 — formula lines mandatory)
    - @packages/HitTriage/src/app/hit-triage-app.ts (lines 384-410 — grok.shell.user.friendlyName precedent for audit)
  </read_first>
  <action>
Create `src/publishing/publish-project.ts`. Standard imports plus all three Wave 1 primitives, `createVolcanoPlot` from `../viewers/volcano`, `FORMULA_LINE_ASSERTION_PREFIX` from `./assert-published-shape`, plus `awaitCheck`/`delay` from `@datagrok-libraries/test/src/test` (note: this is fine even in non-test code per existing package precedent).

**Export `applyVolcanoFormulaLines(viewer: DG.Viewer, fcThreshold: number, pThreshold: number): void`** — pure helper that mutates the passed viewer to add (or replace) the two formula lines. Body:
- Build the formula-line array shape (look at `volcano.ts` for the existing precedent — typically an array of `{title, formula, color, style?}` objects):
  - `{title: 'FC = +' + fcThreshold, formula: '${x} = ' + fcThreshold, color: '#888888'}`
  - `{title: 'FC = -' + fcThreshold, formula: '${x} = -' + fcThreshold, color: '#888888'}`
  - `{title: 'p = ' + pThreshold, formula: '${y} = -log10(' + pThreshold + ')', color: '#888888'}`
- Apply via `viewer.setOptions({look: {formulaLines: JSON.stringify(lines)}})` (canonical platform shape; verify exact key against `volcano.ts` precedent — if volcano uses `props.formulaLines` instead, match that).
- Idempotent — overwrites any existing formula lines.
- Exported because Plan 07's post-open recovery hook re-uses this helper on reopened published Projects (B-2 / W-7 belt-and-braces #3).

Export `publishAnalysis(df: DG.DataFrame, opts: PublishOptions): Promise<DG.Project>`.

Wrap the entire body in `try { ... } finally { pi.close(); }` where `pi = DG.TaskBarProgressIndicator.create('Publishing snapshot...')`. Each step calls `pi.update('<step name>')` so the analyst sees progress per RESEARCH §"Pitfall 14" UX note.

**Step 1 — assemble PublishedMetadata:**
- Read sharer identity SERVER-AUTHORITATIVE per PATTERNS.md Section 3 Analog C: `const publishedBy = grok.shell.user.friendlyName;` and `const publishedByEmail = grok.shell.user.email ?? null;`. NEVER accept a client-supplied `published_by` (Security Mistake #7 per RESEARCH §"Don't Hand-Roll").
- Read DE method + thresholds from source DF tags: `df.getTag('proteomics.de_method')`, plus `fcThreshold` and `pThreshold` (typically stored on the volcano view's options, NOT the DF — read from `grok.shell.tv?.getOptions()` or default to `1.0` / `0.05` if not retrievable; Plan 05 dialog can override these via `opts`).
- Compute slug + version: `const slug = slugifyTarget(opts.target);` and `const version = opts.priorVersion ? <parse prior version + 1> : 1;` (regex `/-v(\d+)-/` on prior's name).
- Generate `publishId = crypto.randomUUID()` (browser-native).
- `dateStr = new Date().toISOString().slice(0, 10);` for the project name.
- Set `meta: PublishedMetadata` with all fields populated. Initially `meta.supersedes = opts.priorVersion?.id ?? null` and `meta.supersededBy = null`.
- Detect enrichment carry: look for an enrichment DF in `grok.shell.tables` via `tables.find(t => t.getTag('proteomics.enrichment') === 'true')`. Set `meta.includesEnrichment = !!enrichSource`.

**Step 2 — trim protein DF:**
- `pi.update('Trimming snapshot...')`
- `const frozen = trimForPublish(df, meta);` (Pitfall 1 mitigation — source `df` is NOT mutated).
- Set `const expectedName = frozen.name;` (Plan 02 sets it to `<source>_published_<date>`).
- Note: `setPublishedTags(frozen, meta)` inside `trimForPublish` writes `PUBLISHED_FC_THRESHOLD` and `PUBLISHED_P_THRESHOLD` tags on `frozen` automatically (Plan 01 task 2 writes every PUBLISHED_TAGS key). This is belt-and-braces #2 — the post-open recovery hook (Plan 07) reads these tags to re-apply formula lines if the serializer strips look config.

**Step 3 — trim enrichment if present (opportunistic, D-05):**
- `let frozenEnrich: DG.DataFrame | null = null;`
- `if (enrichSource) frozenEnrich = trimEnrichmentForPublish(enrichSource, meta);`

**Step 4 — ensure umbrella Space (umbrellaName override support):**
- `pi.update('Ensuring umbrella Space...')`
- `const UMBRELLA_NAME = opts.umbrellaName ?? 'Proteomics-Reviews';` — opts.umbrellaName override supports Plan 08 Test 5 negative test (throwaway umbrella); production callers leave it undefined.
- Per RESEARCH §"Pattern 2" step 3:
  - `if (await grok.dapi.spaces.rootSpaceExists(UMBRELLA_NAME))` then `umbrella = await grok.dapi.projects.filter('name = "' + UMBRELLA_NAME + '" and isSpace = true').first();`
  - else `umbrella = await grok.dapi.spaces.createRootSpace(UMBRELLA_NAME);`
- `const umbrellaClient = grok.dapi.spaces.id(umbrella.id);`

**Step 5 — ensure per-target child Space:**
- `pi.update('Ensuring per-target Space...')`
- `const childName = 'Proteomics-Review-' + slug;`
- Per RESEARCH §"Pattern 2" step 4:
  - `if (await umbrellaClient.subspaceExists(childName))` then list children via `umbrellaClient.children.filter('Project', false).list()` and pick by `friendlyName === childName` (RESEARCH Open Question 3 — no `subspaceById(name)` method).
  - else `childSpace = await umbrellaClient.addSubspace(childName);`
- `const childClient = grok.dapi.spaces.id(childSpace.id);`

**Step 6 — save Project (multi-child pattern per Bio):**
- `pi.update('Saving Project...')`
- `const project = DG.Project.create();`
- `project.name = 'Proteomics-Review-' + slug + '-v' + version + '-' + dateStr;`
- Set `meta.publishId` onto `project.options['proteomics.published_id']` so the orchestrator can re-find this Project on supersede.
- `project.addChild(frozen.getTableInfo());` and if `frozenEnrich`: `project.addChild(frozenEnrich.getTableInfo());`

**Step 6.5 — REVISION B-2: re-create volcano on the trimmed DF + apply formula lines BEFORE addChild(tv.getInfo()) AND BEFORE the save:**
- Open the trimmed DataFrame in its own TableView: `const trimmedTv = grok.shell.addTableView(frozen); await delay(100);`
- Re-create the volcano against the trimmed DF: `const volcano = createVolcanoPlot(frozen);` then `trimmedTv.addViewer(volcano);` (match `src/package.ts` precedent for how volcano is docked; `grok.shell.dockManager.dock(volcano, ...)` if needed).
- Apply formula lines: `applyVolcanoFormulaLines(volcano, meta.fcThreshold, meta.pThreshold);`
- This is belt-and-braces #1 — even if the serializer strips look config on save/reopen, the publish-side persistence path stamped the values into tags (belt-and-braces #2 via setPublishedTags) and the post-open recovery hook (Plan 07) re-applies them on first paint after reopen (belt-and-braces #3 via this same exported `applyVolcanoFormulaLines`).
- Capture `const viewInfo = trimmedTv.getInfo();` for the addChild step below.
- `project.addChild(viewInfo);`
- `await grok.dapi.tables.uploadDataFrame(frozen); await grok.dapi.tables.save(frozen.getTableInfo());`
- If enrichment: same for `frozenEnrich`.
- `await grok.dapi.views.save(viewInfo);`
- `await grok.dapi.projects.save(project);`
- `meta.publishId = project.id;` (update with the server-assigned id, for downstream supersede)
- Move Project into child Space: `await childClient.addEntity(project.id);`

Return the partially-complete `project` at end of this task — Tasks 2, 3, 4 add steps 7, 8, 9. Initial scaffold checks in TypeScript.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/publish-project.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export async function publishAnalysis" src/publishing/publish-project.ts | grep -q '^1$' &amp;&amp; grep -c "export function applyVolcanoFormulaLines" src/publishing/publish-project.ts | grep -q '^1$' &amp;&amp; grep -c "trimForPublish" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "createVolcanoPlot" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "createRootSpace\|rootSpaceExists" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "addSubspace" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "grok\.shell\.user\.friendlyName" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "opts\.umbrellaName" src/publishing/publish-project.ts | grep -qv '^0$'</automated>
  </verify>
  <done>`publishAnalysis` exported; `applyVolcanoFormulaLines` exported as separate helper; reads `grok.shell.user.friendlyName` server-side (never trusts client-supplied); calls `trimForPublish` (and `trimEnrichmentForPublish` opportunistically); ensures umbrella + child Space via D-03 primary shape with `opts.umbrellaName` override support; step 6.5 re-creates volcano + applies formula lines via `applyVolcanoFormulaLines` BEFORE the save; saves multi-child Project via Bio pattern; moves Project into child Space via `addEntity`.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 2: Implement step 7 — permission grant + verify-and-rollback gate (T-15-01)</name>
  <files>src/publishing/publish-project.ts</files>
  <read_first>
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-03 exact error string for verify-and-rollback)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Pattern 2" step 6 + step 7 — verify-and-rollback gate; Assumptions A1, A2)
    - @.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md (spike A1 + A2 outputs — locked permissions.get shape + inheritance behavior)
    - @packages/ApiSamples/scripts/dapi/layouts-and-permissions.js (`permissions.grant(entity, group, edit:bool)` + `permissions.get(entity)` shape)
  </read_first>
  <action>
Add step 7 to `publishAnalysis` after step 6.5 (Task 1):

**Step 7a — grant View to reviewer group at child-Space level (D-03):**
- `pi.update('Granting reviewer access...')`
- `await grok.dapi.permissions.grant(childSpace, opts.reviewerGroup, false);` — the third arg `edit: false` means View-only (per ApiSamples `layouts-and-permissions.js`).
- All Projects already in (and subsequently moved into) the child Space inherit View per Plan 00 spike A2.

**Step 7b — VERIFY-AND-ROLLBACK GATE (Pitfall 2 mitigation, T-15-01, NON-NEGOTIABLE):**
- `pi.update('Verifying view-only access...')`
- `const perm = await grok.dapi.permissions.get(childSpace) as {view: DG.Group[]; edit: DG.Group[]};`
- `const inView = perm.view.some((g) => g.id === opts.reviewerGroup.id);`
- `const inEdit = perm.edit.some((g) => g.id === opts.reviewerGroup.id);`
- If `!inView || inEdit`:
  - **ROLLBACK:** `try { await grok.dapi.projects.delete(project); } catch (rollbackErr) { grok.shell.error('Manual cleanup required: project id ' + project.id + ' could not be deleted: ' + (rollbackErr as any)?.message); }` — T-15-04 mitigation: never let rollback errors silently swallow the original failure.
  - **THROW exact error string per CONTEXT.md D-03:** `throw new Error('Reviewer group already has elevated access via Space inheritance — publish aborted; ask an admin to scope the umbrella Space\\'s permissions');`

This step + step 8 are the TWO non-negotiable gates per CONTEXT.md. There is no "publish succeeded but ACL is wrong" partial-success state.

Per RESEARCH Open Question 6 / Assumption A4: per-target ACL audit may also want `permissions.get(project)` AFTER inheritance — defensive secondary check. If spike A2 confirmed grant inheritance flows childSpace -> contained Project, ONE check on `permissions.get(childSpace)` is sufficient; if NOT, ALSO run `const perm2 = await grok.dapi.permissions.get(project);` and assert reviewer group is in `perm2.view`, NOT in `perm2.edit`. Either way, rollback on any failure.

Add an additional defensive check per CONTEXT.md "non-negotiable": if `Object.keys(perm)` includes anything beyond `view` + `edit` (e.g. `share`, `delete` per RESEARCH A1 note), assert reviewer group is NOT in those either. Rollback if found.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/publish-project.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "permissions\\.grant" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "permissions\\.get" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "projects\\.delete" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -q "Reviewer group already has elevated access via Space inheritance" src/publishing/publish-project.ts &amp;&amp; grep -c "Manual cleanup required" src/publishing/publish-project.ts | grep -qv '^0$'</automated>
  </verify>
  <done>Step 7a grants View (third arg `false`); Step 7b reads `permissions.get(childSpace)`, asserts reviewer in `view` AND not in `edit` (and not in `share`/`delete` if present); on failure deletes project + throws with EXACT error string from D-03; rollback delete is wrapped in try/catch with manual-cleanup surface if it fails (T-15-04).</done>
</task>

<task type="auto" tdd="false">
  <name>Task 3: Implement step 8 — assertPublishedShape gate with self-healing on FORMULA_LINE_MISSING (W-7) + step 9 — supersede chain dual-write (W-8)</name>
  <files>src/publishing/publish-project.ts</files>
  <read_first>
    - @packages/Proteomics/src/publishing/assert-published-shape.ts (just-built helper from Plan 03 + FORMULA_LINE_ASSERTION_PREFIX constant)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-04 supersede chain — bidirectional + prior NOT deleted; PUB-06 SC-2 — formula lines mandatory; "Belt-and-braces is the design philosophy")
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Open Question 6" — both DF tag AND Project.options for supersede pointer)
    - @.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md (spike A4 — Project.options survival on reopen)
  </read_first>
  <action>
Add steps 8 + 9 to `publishAnalysis` after step 7 (Task 2):

**Step 8 — assertPublishedShape ROUND-TRIP GATE with self-healing on FORMULA_LINE_MISSING (Pitfall 3, NON-NEGOTIABLE):**
- `pi.update('Verifying round-trip survival...')`
- Build `PublishedShapeContract` from local state:
  - `expectedName: frozen.name`
  - `expectedProjectName: project.name`
  - `expectedAllowlist: frozen.columns.toList().filter(c => !c.name.startsWith('_meta_')).map(c => c.name)`
  - `expectedMeta: meta`
  - `expectVolcano: true`
  - `expectEnrichment: !!frozenEnrich`
  - `expectFormulaLines: true` (PUB-06 / SC-2 — formula lines are mandatory on the reviewer first paint)

- **Self-healing wrapper around assertPublishedShape (W-7):**
  ```
  let healedOnce = false;
  while (true) {
    try {
      await assertPublishedShape(project, contract);
      break; // success
    } catch (assertErr: any) {
      const msg = String(assertErr?.message ?? '');
      // SELF-HEAL PATH: serializer stripped formula lines — re-apply once via the post-open recovery hook
      if (!healedOnce && msg.includes(FORMULA_LINE_ASSERTION_PREFIX)) {
        healedOnce = true;
        // Re-apply formula lines on the currently-open reopened volcano.
        // assertPublishedShape has just reopened the project into grok.shell.tv; the volcano viewer
        // is in tv.viewers. Find it and re-apply.
        const tv = grok.shell.tv;
        const volcano = tv ? Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT) : null;
        if (volcano) {
          // Read thresholds from the reopened DF tags (belt-and-braces #2 - tag survived even though look config didn't)
          const reDf = tv!.dataFrame;
          const fc = parseFloat(reDf.getTag(PUBLISHED_TAGS.PUBLISHED_FC_THRESHOLD) ?? String(meta.fcThreshold));
          const p = parseFloat(reDf.getTag(PUBLISHED_TAGS.PUBLISHED_P_THRESHOLD) ?? String(meta.pThreshold));
          applyVolcanoFormulaLines(volcano, fc, p);
          // Re-save the view so the recovery sticks for future reopens.
          try { await grok.dapi.views.save(tv!.getInfo()); } catch { /* best-effort */ }
          continue; // retry assertion
        }
      }
      // Either NOT a formula-line failure, OR we already tried to heal once — roll back.
      try { await grok.dapi.projects.delete(project); } catch (rb: any) {
        grok.shell.error('Manual cleanup required: project id ' + project.id + ' could not be deleted: ' + (rb?.message ?? rb));
      }
      throw new Error('Round-trip shape verification failed: ' + (assertErr?.message ?? String(assertErr)));
    }
  }
  ```

After successful assertion (possibly after one self-heal pass), the reviewer-facing artifact is verified-as-shipping.

**Step 9 — supersede chain (D-04 + W-8 dual-write, conditional on republish):**
- `pi.update('Updating supersede chain...')`
- If `opts.priorVersion !== null`:
  - **W-8 belt-and-braces — write `superseded_by` to BOTH the prior Project's options AND the prior Project's DataFrame tag** so Plan 06 panel's "Newer version available" link surfaces on reopen of the prior version regardless of which read path the serializer drops:
    - Prior Project's options: `opts.priorVersion.options['proteomics.superseded_by'] = project.id;` then `await grok.dapi.projects.save(opts.priorVersion);`
    - Prior Project's DataFrame tag — reopen the prior project to get its DF, stamp the tag, re-save:
      ```
      try {
        grok.shell.closeAll();
        await delay(100);
        const priorReopened = await grok.dapi.projects.find(opts.priorVersion.id);
        await priorReopened.open();
        await awaitCheck(() => !!grok.shell.tv?.dataFrame, 'supersede stamp: prior project DF did not materialize', 5000);
        const priorDf = grok.shell.tv!.dataFrame;
        priorDf.setTag(PUBLISHED_TAGS.SUPERSEDED_BY, project.id);
        // belt-and-braces column write if the column already exists (Plan 02 wrote _meta_superseded_by at trim time;
        // it was empty string for v1 because supersededBy was null then; now we patch it).
        const supersededByCol = priorDf.col(META_COLUMNS.SUPERSEDED_BY);
        if (supersededByCol) {
          // Single-row metadata column — overwrite row 0
          supersededByCol.set(0, project.id);
        }
        // Re-save the prior's table info so the tag + column persist.
        await grok.dapi.tables.uploadDataFrame(priorDf);
        await grok.dapi.tables.save(priorDf.getTableInfo());
      } catch (stampErr: any) {
        grok.shell.warning('Supersede tag could not be stamped on prior version DF: ' + (stampErr?.message ?? stampErr) + ' (Project.options pointer was still written — Plan 06 panel falls back to that path)');
      }
      ```
  - On the new Project: `project.options['proteomics.supersedes'] = opts.priorVersion.id;` then `await grok.dapi.projects.save(project);` (note: `frozen` already has `proteomics.supersedes` tag set at trim time via `setPublishedTags(frozen, meta)` when `meta.supersedes = opts.priorVersion.id` was assigned in Step 1).
- Prior version is NEVER deleted (D-04 explicit: "soft pointer, never destructive").

**Return:**
- `pi.update('Done.');` (the finally block closes pi)
- Show success via `grok.shell.info('Published: ' + project.name);`
- `return project;`

**Error wrapping for surface presentation:**
The orchestrator throws specific errors at each step. Plan 05 (`share-dialog.ts`) catches and surfaces via `grok.shell.error(...)`. Add JSDoc note:
- Step 7 throws with the EXACT error string from D-03
- Step 8 throws "Round-trip shape verification failed: <specific cause>" — but only AFTER the self-healing path was tried once for FORMULA_LINE_MISSING failures
- Any earlier step throws on its own (no rollback for steps 1-6 since nothing is server-persisted yet, except the umbrella Space creation in Step 4 which is harmless to leave behind for future publishes)
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/publish-project.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "assertPublishedShape" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "FORMULA_LINE_ASSERTION_PREFIX" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "applyVolcanoFormulaLines" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "healedOnce\\|healed once" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "opts\\.priorVersion" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "superseded_by" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "supersedes" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "Round-trip shape verification failed" src/publishing/publish-project.ts | grep -qv '^0$' &amp;&amp; grep -c "TaskBarProgressIndicator" src/publishing/publish-project.ts | grep -qv '^0$'</automated>
  </verify>
  <done>Step 8 wraps assertPublishedShape in self-healing loop that detects FORMULA_LINE_ASSERTION_PREFIX, re-applies formula lines via applyVolcanoFormulaLines reading thresholds from tags, retries once, and only rolls back if STILL failing (W-7); Step 9 writes bidirectional supersede pointers via DUAL-WRITE (W-8 — Project.options AND DataFrame tag + metadata column on the prior); prior version is never deleted; TaskBarProgressIndicator updates each phase; strict TypeScript compiles.</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| publishing user (analyst) -> reviewer (biologist consumer) | Publishing user produces a frozen artifact; reviewer must only have View, never Edit |
| Datagrok server permissions API -> client orchestrator decision | Client orchestrator trusts `permissions.get` server-truth; if Edit slipped in, rollback |
| live source DF -> trimForPublish clone boundary | Trim primitive owns the clone; orchestrator NEVER mutates source DF directly |
| platform serializer -> reopened look config | Look config (formula lines) may be stripped; orchestrator writes belt-and-braces tags + provides self-heal path + recovery hook |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-01 | Elevation of privilege (Space-inheritance ACL leak, HIGH severity per Pitfall 2) | publishAnalysis step 7 verify-and-rollback gate | mitigate | Non-negotiable post-grant `permissions.get` check; reviewer group MUST be in `view` AND NOT in `edit`/`share`/`delete`; on failure delete Project + throw exact D-03 error string |
| T-15-02 | Information disclosure (stale-snapshot leak, HIGH severity per Pitfall 1) | publishAnalysis step 2 (trim via trimForPublish) | mitigate | `trimForPublish` is the orchestrator's FIRST mutation primitive; source DF never touched by orchestrator (trim creates the clone — see Plan 02 for the per-primitive enforcement). Plan 08 Test 3 asserts source unchanged after publish. (W-9: T-15-02 row explicitly enumerated here at orchestrator level since this orchestrator is the integration point that calls trimForPublish.) |
| T-15-03 | Tampering (freeform target injection into Space name / Project name) | publishAnalysis step 1 (slug) | mitigate | `slugifyTarget` from Plan 01 strips to `[A-Za-z0-9._-]`, caps 64 chars, never empty; raw target preserved in tag/column only |
| T-15-04 | Denial of service (Project deletion on rollback edge case) | publishAnalysis steps 7 + 8 rollback paths | mitigate | Rollback `dapi.projects.delete` wrapped in try/catch; surface "Manual cleanup required: project id X" via `grok.shell.error` if delete fails; never silently swallow the original failure |
| T-15-SP-05 | Tampering (client-supplied `published_by`) | publishAnalysis step 1 metadata assembly | mitigate | Read `grok.shell.user.friendlyName` server-side; reject any caller-supplied `publishedBy` field |
| T-15-SP-12 | Information disclosure (formula lines stripped by serializer — silent loss of SC-2 visual gate) | publishAnalysis steps 6.5 + 8 self-heal | mitigate | Belt-and-braces #1 (apply formula lines pre-save), #2 (persist threshold values in tags via setPublishedTags), #3 (self-heal in step 8 catch + post-open recovery hook in Plan 07) |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- 9-step sequence respected (no step skipped, steps 7 + 8 are non-negotiable gates; step 6.5 applies formula lines pre-save)
- Step 8 self-healing path triggers on FORMULA_LINE_ASSERTION_PREFIX exactly once before rolling back (W-7)
- Step 9 writes superseded_by via dual-write path: Project.options on prior + DataFrame tag on prior (W-8)
- Both rollback paths (step 7 + step 8) delete the Project AND throw with descriptive error
- Manual-cleanup surface for failed rollback (T-15-04)
- `applyVolcanoFormulaLines` exported so Plan 07 post-open recovery hook can re-use it
- `opts.umbrellaName` override supports Plan 08 Test 5 negative test
- Smoke verification (deferred to Plan 08 test suite — round-trip with synthetic fixture)
</verification>

<success_criteria>
- `src/publishing/publish-project.ts` exports `publishAnalysis(df, opts): Promise<DG.Project>` AND `applyVolcanoFormulaLines(viewer, fc, p): void`
- Function never mutates source DF directly (Pitfall 1)
- Function reads `grok.shell.user.friendlyName` server-side (Security Mistake #7 mitigation)
- Step 6.5 re-creates volcano on trimmed DF and applies formula lines BEFORE the save (B-2)
- Step 7 throws EXACT error string from D-03 on inheritance leak
- Step 8 self-heals on FORMULA_LINE_MISSING once, then throws "Round-trip shape verification failed: ..." on persistent failure (W-7)
- Step 9 writes superseded_by to BOTH Project.options AND DataFrame tag (W-8)
- Prior version is NEVER deleted on republish (D-04 explicit)
- TaskBarProgressIndicator updates at each step (Pitfall 14 UX requirement)
- T-15-02 row explicitly enumerated in this plan's STRIDE register (W-9)
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-04-SUMMARY.md` when done with:
- 9-step sequence + step 6.5 confirmed implemented
- Notes on which spike outputs (A1, A2, A4) drove specific code choices
- Confirmation that both rollback gates use `dapi.projects.delete` + descriptive error strings
- Confirmation that supersede chain is bidirectional, dual-write (Project.options + DF tag/column), and non-destructive
- Confirmation that the self-healing path on FORMULA_LINE_MISSING was tested (or noted as a deferred Plan 08 case)
- Confirmation that `applyVolcanoFormulaLines` is exported for Plan 07 post-open recovery hook re-use
- Confirmation that `opts.umbrellaName` override is wired through to step 4
</output>
