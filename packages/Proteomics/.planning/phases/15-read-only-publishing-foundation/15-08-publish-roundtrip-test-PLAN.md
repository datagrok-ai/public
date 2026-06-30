---
phase: 15-read-only-publishing-foundation
plan: 08
type: execute
wave: 6
depends_on: ["15-04", "15-07"]
files_modified:
  - src/tests/publish-roundtrip.ts
autonomous: true
requirements: [PUB-03, PUB-06, PUB-12]
must_haves:
  truths:
    - "`grok test --category Publishing` runs all 5 publish-roundtrip regression tests against the live server and reports green when behavior matches the contract"
    - "Test 1 — assertPublishedShape round-trip on synthetic fixture (no-enrichment shape) — exercises the full publish + reopen + assert cycle AND asserts that the reopened volcano carries formula lines for the FC threshold and p threshold (PUB-06 / B-2)"
    - "Test 2 (P2 — PUB-12) — assertPublishedShape round-trip on Spectronaut Candidates fixture (with-enrichment shape, D-05); asserts that on reopen the enrichment-to-protein cross-DF current-row subscription is re-established and a current-row change in the enrichment DF triggers a protein selection update (W-5)"
    - "Test 3 — mutate source DF after publish; assert reopened snapshot is frozen unchanged (Pitfall 1)"
    - "Test 4 — republish creates a NEW Project with bidirectional `superseded_by`/`supersedes` pointers (D-04); asserts the superseded_by reaches the prior DataFrame's tag as well as Project.options (W-8 dual-write)"
    - "Test 5 — verify-and-rollback gate REJECTS Edit-inherited grant with the EXACT D-03 error string; new Project is deleted (Pitfall 2 negative test)"
    - "Every test cleans up its own server-side entities (projects + spaces) on completion"
  artifacts:
    - path: "src/tests/publish-roundtrip.ts"
      provides: "Load-bearing gate per success criterion 3 — 5 round-trip regression tests"
      contains: "category('Publishing', ...)"
  key_links:
    - from: "src/tests/publish-roundtrip.ts"
      to: "src/publishing/publish-project.ts"
      via: "publishAnalysis call"
      pattern: "publishAnalysis\\("
    - from: "src/tests/publish-roundtrip.ts"
      to: "src/publishing/assert-published-shape.ts"
      via: "assertPublishedShape call"
      pattern: "assertPublishedShape\\("
---

<objective>
Build the load-bearing regression suite that guarantees Phase 15's contract holds across future changes. 5 tests cover the 5 ROADMAP success criteria + the two non-negotiable gates + the B-2 / W-5 belt-and-braces paths added in this revision pass. Anything that breaks any of these tests breaks the phase contract.

Purpose: ROADMAP success criterion 3 names this file explicitly — "the `assertPublishedShape` round-trip test passes." Plan 08 is the criterion. Plan 04's orchestrator calls `assertPublishedShape` defensively in-process; this plan exercises the same helper as a regression gate so a future refactor that subtly breaks the round-trip is caught at `grok test` time, not at production-share time.

**Revision additions:**
- **B-2:** Test 1 explicitly asserts formula-line presence on the reopened volcano (PUB-06 / SC-2). This was missing — covered only by the in-orchestrator `assertPublishedShape` call (which Plan 04 self-heals in process) — the regression suite needs an INDEPENDENT post-`projects.find().open()` check.
- **W-5:** Test 2 explicitly asserts the cross-DF enrichment-to-volcano subscription is re-established on reopen. Drives a current-row change in the enrichment DF, asserts the protein DF's selection updates.
- **Task-level P2:** Test 2 is P2-tagged (PUB-12 enrichment carry is P2 per CONTEXT.md). Tests 1, 3, 4, 5 are P1.

Output: `src/tests/publish-roundtrip.ts` registering 5 tests under `category('Publishing', ...)`. Tests are categorized so `grok test --category Publishing` runs the full suite in ~30-60 seconds.
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@.planning/phases/15-read-only-publishing-foundation/15-VALIDATION.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/tests/analysis.ts
@packages/Proteomics/src/publishing/publish-project.ts
@packages/Proteomics/src/publishing/assert-published-shape.ts
@packages/Proteomics/src/publishing/publish-state.ts
@packages/Proteomics/src/publishing/trim-dataframe.ts
@packages/Proteomics/src/publishing/post-open-recovery.ts
@packages/Bio/src/tests/projects-tests.ts

<interfaces>
<!-- Test framework -->

From @datagrok-libraries/test/src/test:
  category(name: string, fn: () => void): void
  test(name: string, fn: () => Promise<void> | void, options?: {skipReason?: string; timeout?: number}): void
  expect(actual, expected): void
  awaitCheck(predicate, errorMsg, timeoutMs): Promise<void>
  delay(ms): Promise<void>

From Plan 04 (publish-project.ts):
  publishAnalysis(df, opts): Promise<DG.Project>
  applyVolcanoFormulaLines(viewer, fc, p): void

From Plan 03 (assert-published-shape.ts):
  assertPublishedShape(project, contract): Promise<void>
  PublishedShapeContract: { ..., expectFormulaLines: boolean }
  FORMULA_LINE_ASSERTION_PREFIX

From Plan 01 (publish-state.ts):
  PUBLISHED_TAGS, META_COLUMNS, isPublished, getPublishedMetadata, slugifyTarget

From Plan 07 (post-open-recovery.ts):
  recoverPublishedProject(df): Promise<void>
  wireEnrichmentToVolcano(enrichDf, proteinDf): void

From datagrok-api:
  grok.dapi.permissions.grant(entity, group, edit: boolean): Promise<void>
  grok.shell.user.group: DG.Group
  grok.dapi.spaces.createRootSpace(name): Promise<DG.Project>
  grok.dapi.spaces.delete(space): Promise<void>     // cascades to children per spike A3
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Build fixture helpers + Test 1 (no-enrichment round-trip with formula-line assertion — B-2)</name>
  <files>src/tests/publish-roundtrip.ts</files>
  <read_first>
    - @packages/Proteomics/src/tests/analysis.ts (lines 1-60 — in-package test convention)
    - @packages/Bio/src/tests/projects-tests.ts (lines 26-97 — save/reopen + assertion shape)
    - @packages/Proteomics/src/publishing/publish-project.ts (publishAnalysis signature + behavior + applyVolcanoFormulaLines)
    - @packages/Proteomics/src/publishing/assert-published-shape.ts (PublishedShapeContract structure including expectFormulaLines)
    - @packages/Proteomics/src/publishing/trim-dataframe.ts (trim contract for expected fixture column list)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 7 — test analog + 7 required test cases)
    - @.planning/phases/15-read-only-publishing-foundation/15-VALIDATION.md (load-bearing test enumeration)
  </read_first>
  <action>
Create `src/tests/publish-roundtrip.ts`. Standard imports + `category`, `test`, `expect`, `awaitCheck`, `delay` from `@datagrok-libraries/test/src/test` + all Phase 15 surfaces.

**Fixture helpers:**

`createSyntheticDeFixture(): DG.DataFrame` — builds a 10-row no-enrichment fixture. Match the schema Plan 02 trim expects:
- Columns: `Protein ID` (string), `Gene Name` (string), `log2FC` (float), `p-value` (float), `adj.p-value` (float), `significant` (string `'yes'`/`'no'`), `direction` (string `'Up'`/`'Down'`/`'NS'`).
- 10 deterministic rows (5 significant Up, 3 significant Down, 2 NS). Compute `adj.p-value` from `p-value * 10 / rank` style — anything reproducible.
- Tags: `proteomics.source='generic'`, `proteomics.de_method='limma'`, `proteomics.de_complete='true'`, `proteomics.groups=JSON.stringify({group1:{name:'Control',columns:['s1','s2','s3']},group2:{name:'Treatment',columns:['s4','s5','s6']}})`.
- SemTypes assigned on `Protein ID`, `Gene Name`, `log2FC`, etc.
- `df.name = 'spike-fixture-no-enrichment'`.

`createSpectronautCandidatesFixture(): {protein: DG.DataFrame, enrichment: DG.DataFrame}` — builds the with-enrichment shape:
- protein DF: same shape as above but `proteomics.source='spectronaut-candidates'`, `proteomics.de_method='spectronaut'`, `proteomics.enrichment='true'` (marker for D-05).
- enrichment DF: 6 rows with `Term Name`, `Source` (`'GO:BP'`, `'KEGG'`), `p-value`, `adj.p-value`, `Intersection` (the gene-list column wireEnrichmentToVolcano subscribes to — populate with gene names that exist in the protein DF so the subscription test in Test 2 can verify the selection update). Tag: `proteomics.enrichment='true'`.

`pickAnyGroup(): Promise<DG.Group>` — `return grok.shell.user.group;` (the user's own group is always available + a safe test subject; per Plan 04 only the verify-and-rollback gate cares about WHICH group).

`buildExpectedContract(frozen: DG.DataFrame, project: DG.Project, meta: PublishedMetadata, expectEnrichment: boolean): PublishedShapeContract` — pure helper assembling the contract from local state. Sets `expectFormulaLines: true` by default (PUB-06 / SC-2 mandate).

`cleanupProject(project: DG.Project): Promise<void>` — wraps `try { await grok.dapi.projects.delete(project); } catch { /* swallow — best-effort */ }`.

**`category('Publishing', () => { ... })` registration.**

**Test 1 — `'assertPublishedShape round-trip on synthetic fixture (no-enrichment) — formula lines survive reopen'`** (P1; covers PUB-01, PUB-02, PUB-03, PUB-04, PUB-06, PUB-07, PUB-09, PUB-11):
- Build synthetic fixture via `createSyntheticDeFixture()`.
- Open in TableView: `const tv = grok.shell.addTableView(df); await delay(100);`
- Add a scatter plot (the simulated volcano): `tv.scatterPlot({x: 'log2FC', y: 'adj.p-value'}); await delay(100);`
- Call `const opts: PublishOptions = { target: 'MYH7-test-' + Date.now(), reviewerGroup: await pickAnyGroup(), note: 'Phase 15 roundtrip', priorVersion: null };`
- `const project = await publishAnalysis(df, opts);`
- `expect(!!project, true);` — publishAnalysis returns the project on success.
- `assertPublishedShape` has already been called INSIDE publishAnalysis (Plan 04 step 8) — so the publish succeeding IS the assertion success. Re-run it independently to surface any post-step-8 regressions:
  - Re-build contract via `buildExpectedContract` with `expectFormulaLines: true`.
  - `await assertPublishedShape(project, contract);` — should not throw.

**B-2 explicit formula-line assertion (independent of assertPublishedShape):**
- After the assertion, the project is open in `grok.shell.tv`. Find the volcano viewer: `const volcano = Array.from(grok.shell.tv!.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT);`
- `expect(!!volcano, true);`
- Read formula lines: `const opts = (volcano as any).getOptions?.(); const fl = opts?.look?.formulaLines ?? (volcano as any).props?.formulaLines;`
- Parse and assert at least one formula references the FC threshold and one references the p threshold:
  ```
  const lines = fl ? (typeof fl === 'string' ? JSON.parse(fl) : fl) : [];
  const fcThreshold = 1.0; // the value publishAnalysis would have stored (default from tv.getOptions or fallback)
  const pThreshold = 0.05;
  const hasFcLine = lines.some((l: any) => String(l.formula ?? '').includes(String(fcThreshold)));
  const hasPLine = lines.some((l: any) => String(l.formula ?? '').includes(String(pThreshold)) || String(l.formula ?? '').includes('log10'));
  expect(hasFcLine, true);
  expect(hasPLine, true);
  ```
- Cleanup: `await cleanupProject(project);`
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/tests/publish-roundtrip.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "category('Publishing'" src/tests/publish-roundtrip.ts | grep -q '^1$' &amp;&amp; grep -c "createSyntheticDeFixture" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "createSpectronautCandidatesFixture" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "publishAnalysis" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "assertPublishedShape" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "formulaLines\\|formula lines" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -cE "test\\('assertPublishedShape round-trip" src/tests/publish-roundtrip.ts | grep -qv '^0$'</automated>
  </verify>
  <done>File exists; fixture helpers + Test 1 registered under `category('Publishing', ...)`; Test 1 explicitly asserts formula-line presence on the reopened volcano (B-2); each test cleans up its own project; strict TypeScript compiles.</done>
</task>

<task type="auto" priority="P2" tdd="false">
  <name>Task 2 (P2 — PUB-12): Test 2 (with-enrichment round-trip + cross-DF subscription re-established assertion — W-5)</name>
  <files>src/tests/publish-roundtrip.ts</files>
  <read_first>
    - @packages/Proteomics/src/publishing/post-open-recovery.ts (wireEnrichmentToVolcano signature + behavior)
    - @packages/Proteomics/src/viewers/enrichment-viewers.ts (current-row subscription pattern)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-05 enrichment carry — opportunistic)
  </read_first>
  <action>
**Task-level P2:** PUB-12 (enrichment carry) is P2 per CONTEXT.md, so this test is P2. Tests 1, 3, 4, 5 are P1.

Add inside the same `category('Publishing', () => { ... })`:

**Test 2 — `'assertPublishedShape round-trip on Spectronaut Candidates fixture (with-enrichment); cross-DF subscription re-established on reopen'`** (P2; covers PUB-12 + D-05 + W-5):

- Build `{protein, enrichment}` via `createSpectronautCandidatesFixture()`.
- Open BOTH in shell: `grok.shell.addTableView(protein); grok.shell.addTableView(enrichment); await delay(100);`
- Make protein the active view (publishAnalysis uses `grok.shell.tv`); add scatter plot.
- Run `publishAnalysis(protein, opts)` with `target: 'CK-test-' + Date.now()`.
- `assertPublishedShape(project, contract)` with `expectEnrichment: true`, `expectFormulaLines: true`.
- Verify the published Project contains 2 DataFrames: query `grok.shell.tables` after reopen, expect at least one with `proteomics.enrichment === 'true'`.

**W-5 cross-DF subscription assertion (after reopen):**
- After `assertPublishedShape` returns (which closed all + reopened), the project is open and the autostart hook (Plan 07) should have re-wired the subscription via `wireEnrichmentToVolcano`. Verify:
  - Find the reopened protein DF and enrichment DF: `const reProteinDf = grok.shell.tables.find((t) => t.getTag(PUBLISHED_TAGS.PUBLISHED) === 'true' && t.getTag('proteomics.enrichment') !== 'true');` and the matching enrichment DF.
  - `expect(!!reProteinDf, true); expect(!!reEnrichDf, true);`
  - Allow the autostart hook ~300ms to fire after view materialization: `await delay(300);`
  - Sentinel check (Plan 07 sets `proteomics.enrichment_wired='true'` after wiring): `expect(reEnrichDf!.getTag('proteomics.enrichment_wired'), 'true');`
  - Drive a current-row change: capture `reProteinDf.selection.trueCount` BEFORE; `reEnrichDf.currentRowIdx = 0; await delay(100);`; capture `reProteinDf.selection.trueCount` AFTER.
  - Assert AFTER > BEFORE (the subscription fired and selected at least one protein matching the enrichment row's gene list).
  - If the assertion fails: also manually call `wireEnrichmentToVolcano(reEnrichDf!, reProteinDf!)` as a fallback to verify the helper itself is correct in isolation; if THAT works, the failure is the autostart hook not firing on time — log a warning and retry once.

- Cleanup: `await cleanupProject(project);`
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/tests/publish-roundtrip.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "with-enrichment\\|Spectronaut Candidates fixture" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "wireEnrichmentToVolcano\\|enrichment_wired" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "currentRowIdx" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "selection\\.trueCount\\|selection\\.init" src/tests/publish-roundtrip.ts | grep -qv '^0$'</automated>
  </verify>
  <done>Test 2 registered as P2 (PUB-12); asserts enrichment DF presence + that the cross-DF subscription is wired (sentinel tag OR drive currentRowIdx and assert protein selection update); strict TypeScript compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 3: Test 3 (Pitfall 1 source-mutation isolation) + Test 4 (republish bidirectional supersede with dual-write — W-8)</name>
  <files>src/tests/publish-roundtrip.ts</files>
  <read_first>
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-04 republish chain — bidirectional + prior NOT deleted)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Pitfall 1" — stale-snapshot leak)
    - @packages/Proteomics/src/publishing/publish-project.ts (Plan 04 step 9 supersede chain — where pointers land: Project.options + DF tag + metadata column on prior — W-8 dual-write)
  </read_first>
  <action>
Add inside the same `category('Publishing', () => { ... })`:

**Test 3 — `'source mutation after publish leaves clone unchanged (Pitfall 1)'`** (P1; PUB-03):
- Build synthetic fixture; publish it; capture `project1.id` and the `frozen` DF's row count (PRE-mutation baseline).
- Reopen project to get the published clone: `const reopened1 = await grok.dapi.projects.find(project1.id); await reopened1.open();`
- Read the clone's column values for `log2FC` of row 0.
- Close all views; switch back to the source DF; mutate it:
  - `df.col('log2FC')!.set(0, 999.999);` (deliberate corruption)
  - `df.columns.remove('Gene Name');` (drop a column)
  - `df.setTag(PUBLISHED_TAGS.PUBLISHED_TARGET, 'CHANGED-TARGET');` (overwrite a tag)
- Reopen the project AGAIN: `const reopened2 = await grok.dapi.projects.find(project1.id); await reopened2.open();`
- Assert clone is UNCHANGED:
  - `expect(grok.shell.tv!.dataFrame.col('log2FC')!.get(0), <original baseline>);` — NOT 999.999
  - `expect(!!grok.shell.tv!.dataFrame.col('Gene Name'), true);` — column still present
  - `expect(grok.shell.tv!.dataFrame.getTag(PUBLISHED_TAGS.PUBLISHED_TARGET), <original target>);` — NOT 'CHANGED-TARGET'
- Cleanup: `await cleanupProject(project1);`

**Test 4 — `'republish creates v2 with bidirectional superseded_by + supersedes pointers (W-8 dual-write)'`** (P1; PUB-10):
- Build synthetic fixture.
- First publish: `const projectV1 = await publishAnalysis(df, opts);` with `opts.priorVersion = null`. Capture `v1Id = projectV1.id; v1Name = projectV1.name;`. Verify `v1Name` matches `/Proteomics-Review-.*-v1-\d{4}-\d{2}-\d{2}/`.
- Second publish (republish): `const opts2 = { ...opts, priorVersion: projectV1 };` then `const projectV2 = await publishAnalysis(df, opts2);`.
- Verify name: `projectV2.name` matches `.*-v2-.*`.
- Reopen `projectV2` and read `project.options['proteomics.supersedes']` — expect `v1Id`.
- Reopen `projectV1`:
  - read `project.options['proteomics.superseded_by']` — expect `projectV2.id` (W-8 Project.options write path).
  - **W-8 dual-write assertion:** ALSO read the reopened prior DF's tag and metadata column:
    - `expect(grok.shell.tv!.dataFrame.getTag(PUBLISHED_TAGS.SUPERSEDED_BY), projectV2.id);` (W-8 DF tag write path)
    - `expect(grok.shell.tv!.dataFrame.col(META_COLUMNS.SUPERSEDED_BY)?.get(0), projectV2.id);` (W-8 metadata column write path — column was empty string at v1 publish, patched on republish per Plan 04 step 9)
- Verify prior version was NOT deleted: `expect(await grok.dapi.projects.find(v1Id), <non-null>);` (resolves to the project entity).
- Cleanup: delete both projects.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/tests/publish-roundtrip.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "source mutation after publish leaves clone unchanged" src/tests/publish-roundtrip.ts | grep -q '^1$' &amp;&amp; grep -c "republish creates v2" src/tests/publish-roundtrip.ts | grep -q '^1$' &amp;&amp; grep -c "supersedes" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "superseded_by\\|SUPERSEDED_BY" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "META_COLUMNS\\." src/tests/publish-roundtrip.ts | grep -qv '^0$'</automated>
  </verify>
  <done>2 additional tests: Pitfall 1 source-mutation isolation (asserts clone unchanged after source mutation) + republish bidirectional supersede (asserts v2 created + Project.options + DF tag + metadata column all carry the pointer per W-8 dual-write + prior NOT deleted); strict TypeScript compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 4: Test 5 (verify-and-rollback rejects Edit-inherited grant — negative test)</name>
  <files>src/tests/publish-roundtrip.ts</files>
  <read_first>
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (D-03 EXACT error string: "Reviewer group already has elevated access via Space inheritance — publish aborted; ask an admin to scope the umbrella Space's permissions")
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Pattern 2" step 7 — verify-and-rollback gate; Assumption A2 inheritance behavior)
    - @packages/Proteomics/src/publishing/publish-project.ts (Plan 04 step 7 — exact throw + rollback; opts.umbrellaName override support)
  </read_first>
  <action>
Add inside the same `category('Publishing', () => { ... })`:

**Test 5 — `'verify-and-rollback rejects Edit-inherited grant with exact D-03 error string'`** (P1; PUB-05, Pitfall 2):

This is a negative test that PROVES the verify-and-rollback gate fires. Setup is: pre-create a throwaway umbrella Space, grant the test group EDIT permission on the umbrella, then attempt `publishAnalysis` with `opts.umbrellaName = <throwaway-name>` (Plan 04 supports this override). The orchestrator's step 4 finds (not creates) the throwaway umbrella; step 5 creates the child Space; step 7 reads `permissions.get(childSpace)` and sees the group has `edit` inherited from umbrella. Gate fires.

Implementation:
- Create throwaway umbrella: `const testUmbrella = await grok.dapi.spaces.createRootSpace('Test-Proteomics-Reviews-' + Date.now());`
- Grant Edit to current user's group ON THE UMBRELLA: `await grok.dapi.permissions.grant(testUmbrella, grok.shell.user.group, true);` (third arg `true` = edit).
- Build fixture; attempt publish with `opts.umbrellaName = testUmbrella.name` (Plan 04 supports this override per revision).
- Assert the publish THROWS with the exact D-03 error string:
  ```
  let caughtError: Error | null = null;
  try {
    await publishAnalysis(df, {...opts, umbrellaName: testUmbrella.name});
  } catch (e: any) {
    caughtError = e;
  }
  expect(caughtError !== null, true);
  expect(caughtError!.message.includes('Reviewer group already has elevated access via Space inheritance'), true);
  ```
- Verify NO new Project was left behind: `grok.dapi.projects.filter('name like "Proteomics-Review-...-v1-..."').list()` — expect 0 matches for the test target.
- Cleanup: delete the throwaway umbrella (cascade per RESEARCH §"Nested Spaces Verification" point 4): `await grok.dapi.spaces.delete(testUmbrella);`

End of category block. Final structural check: 5 tests registered, all under `category('Publishing', ...)`, each with explicit cleanup. Plan 04 umbrellaName-override is no longer a deferred dependency — it's part of the Plan 04 revision and Test 5 runs unskipped.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/tests/publish-roundtrip.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "verify-and-rollback rejects Edit-inherited grant" src/tests/publish-roundtrip.ts | grep -q '^1$' &amp;&amp; grep -c "Reviewer group already has elevated access via Space inheritance" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; grep -c "umbrellaName" src/tests/publish-roundtrip.ts | grep -qv '^0$' &amp;&amp; TEST_COUNT=$(grep -cE "^\\s+test\\(" src/tests/publish-roundtrip.ts) &amp;&amp; [ "$TEST_COUNT" -ge 5 ]</automated>
  </verify>
  <done>5 tests registered under category('Publishing'); Test 5 asserts the exact D-03 error string using `opts.umbrellaName` override (Plan 04 revision); strict TypeScript compiles.</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| test fixture -> live Datagrok server | Tests create + delete real entities; cleanup is the trust contract |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-01 | Elevation of privilege (Pitfall 2 negative test ensures gate fires) | Test 5 verify-and-rollback | mitigate | Test 5 EXPLICITLY verifies the gate fires AND that the new Project is deleted (no partial-success leak); this is the regression guarantee for T-15-01 |
| T-15-02 | Information disclosure (Pitfall 1 negative test ensures clone independence) | Test 3 source mutation isolation | mitigate | Test 3 EXPLICITLY mutates the source DF and verifies the clone is unchanged on reopen; regression guarantee for T-15-02 |
| T-15-04 | Denial of service (leaked test entities) | All tests | mitigate | Every test ends with `cleanupProject` (`try { ... delete } catch { swallow }`); Test 5 ALSO deletes its throwaway umbrella (cascade per spike A3); operator audits via `grok s projects list` post-suite |
| T-15-SP-12 | Information disclosure (formula lines absent — Test 1 explicit assert) | Test 1 formula-line assertion | mitigate | Test 1 EXPLICITLY reads the reopened volcano's formula-line config and asserts FC + p threshold values are present; regression guarantee for B-2 / PUB-06 SC-2 |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- 5 tests registered under category('Publishing')
- Test 1 asserts formula-line presence on the reopened volcano (B-2)
- Test 2 asserts the cross-DF enrichment subscription is re-established on reopen (W-5)
- Test 4 asserts W-8 dual-write (Project.options + DF tag + metadata column on prior)
- Test 5 uses `opts.umbrellaName` override (Plan 04 revision)
- `grok test --category Publishing` runs all 5 against the live server (after `npm run build`)
- All tests clean up their entities (no leaked projects in `grok s projects list` post-run)
</verification>

<success_criteria>
- `src/tests/publish-roundtrip.ts` exists with 5 tests
- Each test corresponds to a row in VALIDATION.md "Load-bearing tests"
- Test 1 verifies formula-line presence on reopen (PUB-06 / SC-2 / B-2)
- Test 2 verifies enrichment subscription is re-wired on reopen (PUB-12 / W-5)
- Test 4 verifies supersede dual-write (PUB-10 / W-8)
- Test 5 exercises the EXACT D-03 error string via `opts.umbrellaName` override
- Suite passes against the live server (deferred to verify-phase)
- Task 2 is P2-tagged at task level (PUB-12 is P2)
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-08-SUMMARY.md` when done with:
- 5 tests + which requirement each covers + which revision finding (B-2, W-5, W-8) each new assertion addresses
- Confirmation Test 5 runs unskipped (Plan 04 umbrellaName override landed in revision)
- Confirmation that `grok test --category Publishing` resolves and runs end-to-end
- Note on cleanup hygiene: no leftover entities in server post-suite
- Task 2 P2 tag confirmed at task level (per B-3 task-level tagging for PUB-12)
</output>
