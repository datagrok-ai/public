---
phase: 15-read-only-publishing-foundation
plan: 00
type: execute
wave: 0
depends_on: []
files_modified:
  - src/tests/publish-spike.ts
  - src/package-test.ts
autonomous: false
requirements: [PUB-02, PUB-06, PUB-11]
must_haves:
  truths:
    - "Running `grok test --category Publishing-Spike` produces a printed enumeration of which `proteomics.*` tags, `Proteomics-*` semTypes, `df.name`, `Project.options` entries, and viewer dock positions survive `DG.Project.save` -> `find(id).open()` against the live server"
    - "Operator can read the spike output and lock the trim allowlist + belt-and-braces metadata column list before Wave 1 commits the trim contract"
    - "Spike runs ONCE, prints shapes, then cleans up its own test project — no leftover entities in the live server"
  artifacts:
    - path: "src/tests/publish-spike.ts"
      provides: "One-shot enumeration probe (Wave 0 only)"
      contains: "category('Publishing-Spike', ...)"
    - path: "src/package-test.ts"
      provides: "Spike test registration so `grok test` picks it up"
      contains: "import './tests/publish-spike';"
  key_links:
    - from: "src/tests/publish-spike.ts"
      to: "live Datagrok server"
      via: "grok.dapi.projects.save + find + open + delete"
      pattern: "DG.Project.create.*projects.save.*projects.find.*open"
---

<objective>
Resolve RESEARCH.md assumptions A1, A2, A3, A4, A6, A8 against the live Datagrok server (`release/1.27.3`) by running a one-shot probe that publishes a minimal fixture, reopens it from id, and prints the actual round-trip survival shape for tags / semTypes / `df.name` / `Project.options` / viewer dock positions.

Purpose: lock the trim contract + belt-and-braces metadata column list BEFORE Wave 1 commits to specific columns. Phase 13 commit `e527d07ba1` already proved the platform serializer partially strips viewer-config tags; this spike enumerates exactly what survives for our specific shape so the Wave 1 trim helper is implemented against verified ground truth, not assumption.

Output: console-printed enumeration that the operator reads once and uses to inform Wave 1 plan execution. The spike does NOT ship as a regression test — it lands in category `Publishing-Spike` (excluded from CI after Wave 0) and self-cleans the test project on completion.
</objective>

<execution_context>
@$HOME/.claude/get-shit-done/workflows/execute-plan.md
@$HOME/.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/PROJECT.md
@.planning/ROADMAP.md
@.planning/STATE.md
@.planning/REQUIREMENTS.md
@.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md
@.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md
@.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md
@.planning/phases/15-read-only-publishing-foundation/15-VALIDATION.md
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/package-test.ts
@packages/Proteomics/src/tests/volcano.ts
@packages/Proteomics/src/utils/proteomics-types.ts

<interfaces>
<!-- Pulled from packages/Bio/src/tests/projects-tests.ts:26-47 — canonical save→find→open round-trip -->
<!-- Spike replaces assertions with console.log of every survivable shape -->

From @datagrok-libraries/test/src/test:
  category(name: string, fn: () => void): void
  test(name: string, fn: () => Promise<void> | void): void
  awaitCheck(predicate: () => boolean, errorMessage: string, timeoutMs: number): Promise<void>
  delay(ms: number): Promise<void>

From datagrok-api:
  DG.Project.create(): DG.Project
  project.options: MapProxy (key-value store on Project, distinct from DataFrame.tags)
  project.addChild(child: DG.Entity): void
  project.name: string
  grok.dapi.tables.uploadDataFrame(df: DG.DataFrame): Promise<void>
  grok.dapi.tables.save(tableInfo: DG.TableInfo): Promise<DG.TableInfo>
  grok.dapi.views.save(layoutInfo: DG.ViewInfo): Promise<DG.ViewInfo>
  grok.dapi.projects.save(project: DG.Project): Promise<DG.Project>
  grok.dapi.projects.find(id: string): Promise<DG.Project>
  grok.dapi.projects.delete(project: DG.Project): Promise<void>
  grok.shell.tv.dataFrame: DG.DataFrame  // active TableView's DataFrame
  grok.shell.tv.viewers: Iterable<DG.Viewer>
  grok.shell.closeAll(): void
</interfaces>
</context>

<tasks>

<task type="auto">
  <name>Task 1: Create publish-spike.ts probe</name>
  <files>src/tests/publish-spike.ts</files>
  <read_first>
    - @packages/Proteomics/src/tests/volcano.ts (in-package test convention — category/test shape)
    - @packages/Bio/src/tests/projects-tests.ts (lines 26-47 saveAndOpenProject + lines 75-97 runSaveAndOpenProjectTest — canonical multi-child save/reopen)
    - @packages/Proteomics/src/utils/proteomics-types.ts (SEMTYPE constants)
    - @packages/Proteomics/src/package.ts (lines 196-203 for volcano viewer type detection precedent)
    - @packages/Proteomics/CLAUDE.md (function-naming + tag conventions)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Pitfall 3 Round-Trip Enumeration" spike outline + Assumptions Log A1, A2, A3, A4, A6, A8)
  </read_first>
  <action>
Create file `src/tests/publish-spike.ts` that registers a single test under `category('Publishing-Spike', ...)`. The test must:

(1) Build a minimal fixture DataFrame in-memory (10 proteins, columns: `Protein ID` (string), `Gene Name` (string), `log2FC` (float), `p-value` (float), `adj.p-value` (float), `significant` (string), `direction` (string)). Tag it with the full proteomics namespace we expect to ship in Wave 1: `proteomics.source='generic'`, `proteomics.de_method='limma'`, `proteomics.de_complete='true'`, `proteomics.groups=JSON.stringify({group1:{name:'Control',columns:[]},group2:{name:'Treatment',columns:[]}})`, `proteomics.published='true'`, `proteomics.published_at=<ISO>`, `proteomics.published_by=<friendly>`, `proteomics.published_target='SPIKE-TARGET-A'`, `proteomics.published_de_method='limma'`, `proteomics.published_fc_threshold='1.0'`, `proteomics.published_p_threshold='0.05'`, `proteomics.published_version='1'`, `proteomics.published_id=<uuid>`, `proteomics.published_includes_enrichment='false'`. Assign `Protein ID` column `semType = SEMTYPE.PROTEIN_ID`, `Gene Name` column `semType = SEMTYPE.GENE_SYMBOL`. Set `df.name = 'spike-source'`.

(2) Open it in a TableView (`grok.shell.addTableView(df)`) and add a scatter plot viewer mimicking a volcano (`tv.scatterPlot({x:'log2FC', y:'adj.p-value', size:10})`). Dock it at `FilterPosition.Right` (or whatever the volcano docks at — match `src/package.ts` volcano-add precedent).

(3) Create a `DG.Project`, name it `spike-publish-roundtrip-<timestamp>`, set TWO `project.options` entries (one named `proteomics.superseded_by='dummy-id-A'`, one named `custom_option='hello'`) to test A4. `addChild(tableInfo)`, `addChild(layoutInfo = tv.getInfo())`. Then sequence `uploadDataFrame -> tables.save -> views.save -> projects.save` per the Bio analog.

(4) Capture `project.id`, call `grok.shell.closeAll()`, then `await grok.dapi.projects.find(projId)` and `.open()`.

(5) After reopen, ENUMERATE and `console.log` (use `grok.shell.info` AND `console.log` both, since the test harness sometimes swallows console.log):
  - `df.name` (Assumption A5 — expected SURVIVE)
  - For every `proteomics.*` tag we set: `reDf.getTag(...)` value (resolves Pitfall 3 unknown for our tag namespace)
  - For each column we tagged: `reDf.col(<name>)?.semType` (resolves whether detectors re-fire on reopen)
  - `reopenedProject.options['proteomics.superseded_by']` and `reopenedProject.options['custom_option']` (Assumption A4)
  - `grok.shell.tv.viewers` mapped to `(v) => ({type: v.type, props: v.getOptions?.()})` — captures whether the scatter plot survived AND whether its `x`/`y`/`color` bindings survived (Assumption A6)
  - `grok.dapi.permissions.get(reopenedProject)` shape: print `Object.keys(perm)` AND `JSON.stringify(perm)` (Assumption A1 — resolves whether perm response has fields beyond `view` and `edit`)

(6) Probe Assumption A2 (Space-inheritance): create a throwaway umbrella Space `Spike-Reviews-<timestamp>`, an inheritance child Space `Spike-Review-Inherit-<timestamp>`, save a second copy of the project into the child Space via `spaces.id(childSpace.id).addEntity(reopenedProject.id)`, grant View on the umbrella Space (NOT the child) to a known group (use the publishing user's own group as a safe test subject), then print `permissions.get(childProject)` to surface whether the inherited grant is visible. Cleanup: delete the child Space — `delete` cascades per RESEARCH §"Nested Spaces Verification" point 4.

(7) Probe Assumption A8 (smart-filter `like` syntax): immediately after `projects.save`, run `grok.dapi.projects.filter('name like "spike-publish-roundtrip-%"').list()` and print result count + first match name. Also try `grok.dapi.projects.filter('name contains "spike-publish-roundtrip"').list()` for comparison.

(8) Probe Assumption A3 (delete cascade): capture `tableInfo.id` BEFORE delete; after `dapi.projects.delete(reopenedProject)`, try `await grok.dapi.tables.find(tableInfoId)` and print whether it resolves or throws 404 / returns null. Wrap in try/catch and log the error shape.

(9) Cleanup at end: `await grok.dapi.projects.delete(reopenedProject)` (idempotent — ok if A3 cascade already wiped it; catch errors).

Test name: `'save+open roundtrip — enumerate survivors'`. Use `await delay(100)` between `addTableView` and viewer-add, between operations involving the shell, per the Bio convention. Use `await awaitCheck(...)` after `reopened.open()` to wait for `grok.shell.tv?.dataFrame` to be defined before reading from it (5000ms timeout per Bio precedent).

NO assertions — this is a probe, not a regression test. The output is read by the operator and informs Wave 1 trim contract.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; grok test --category Publishing-Spike --test "save+open roundtrip — enumerate survivors" 2>&amp;1 | tee /tmp/spike-output.log &amp;&amp; grep -E "(df.name survived|de_method tag|PROTEIN_ID semType|options.custom_option|viewers:|Object.keys\(perm\)|filter list result|delete cascade)" /tmp/spike-output.log</automated>
  </verify>
  <done>File `src/tests/publish-spike.ts` exists, defines `category('Publishing-Spike', ...)`, registers exactly one test `save+open roundtrip — enumerate survivors`, runs end-to-end without throwing, and emits printed output covering all 8 enumeration points above. Operator reads the output and is unblocked to proceed to Wave 1.</done>
</task>

<task type="auto">
  <name>Task 2: Register spike test in package-test.ts</name>
  <files>src/package-test.ts</files>
  <read_first>
    - @packages/Proteomics/src/package-test.ts (existing import block, lines 4-20)
    - @packages/Proteomics/.claude/feedback_grok_test_skipbuild_stale.md (memory: rebuild after adding new test files)
  </read_first>
  <action>
Append `import './tests/publish-spike';` at the end of the import block (after the last existing `import './tests/...';` line). Do NOT add `import './tests/publish-roundtrip';` here — that file does not exist yet (Plan 08 creates it). Match existing single-quote + semicolon convention.

After adding the import, the spike test is reachable via `grok test --category Publishing-Spike` once the bundle is rebuilt. Note for executor: per saved memory `feedback_grok_test_skipbuild_stale.md`, `--skip-build` reuses stale test bundles — run `npm run build` (NOT `grok test --skip-build`) before invoking the spike test for the first time.
  </action>
  <verify>
    <automated>grep -c "import './tests/publish-spike';" packages/Proteomics/src/package-test.ts | grep -q '^1$'</automated>
  </verify>
  <done>`src/package-test.ts` contains exactly one `import './tests/publish-spike';` line; the spike test is registered and discoverable by `grok test --category Publishing-Spike`.</done>
</task>

<task type="checkpoint:human-verify" gate="blocking">
  <name>Task 3: Operator reads spike output and confirms Wave 1 trim contract</name>
  <what-built>The probe in Task 1 + Task 2 emits printed output to the test harness. The operator (you) reads it once and records the resolved values for Assumptions A1, A2, A3, A4, A6, A8 in this phase's commit message or in a brief note inside the phase directory.</what-built>
  <how-to-verify>
1. Rebuild the package (`npm run build` from `packages/Proteomics/`) — NEVER use `grok test --skip-build` for the first run after adding a new test file (memory: stale-bundle gotcha).
2. Run `grok test --category Publishing-Spike` against the running Datagrok server.
3. Read the printed output. For each of the 8 enumeration points in Task 1, capture which shape was observed:
   - Which `proteomics.*` tags SURVIVED reopen? (Wave 1 trim helper must explicitly re-set any that did NOT survive even though we set them on the source.)
   - Which `Proteomics-*` semTypes SURVIVED reopen? (If `PROTEIN_ID` semType is null on reopen, detectors did NOT re-fire — Wave 1 may need to re-assign semTypes post-clone.)
   - `project.options['proteomics.superseded_by']` value on reopen? (Resolves whether Plan 04 uses `project.options` OR `df.setTag` for the supersede pointer — RESEARCH Open Question 6.)
   - `permissions.get(...)` shape: what are the keys? Just `view` + `edit`, or also `share` + `delete`? (Locks the verify-and-rollback assertion in Plan 04.)
   - `like` smart-filter result count vs `contains` result count? (Locks `findPriorShare` filter syntax in Plan 01.)
   - Did the volcano viewer survive? Did its `x`, `y`, `color` props survive? (If no — Plan 03 + Plan 04 add a post-open re-init step.)
   - Did `projects.delete` cascade to `tables.find(tableInfoId)`? (If no — Plan 04 may need explicit `tables.delete(...)` in rollback path.)
   - For the inheritance probe: does `permissions.get(childProject)` show the umbrella-Space's grant? (Locks the verify-and-rollback gate's exact assertion in Plan 04.)
4. Type `approved` to proceed, OR type the observed shape (e.g. "tags all survived; semTypes lost; project.options survived; permissions.get returns view+edit only; like works; volcano survived; cascade works; inheritance visible") so the executor of Plan 01+ knows the locked contract.
  </how-to-verify>
  <resume-signal>Type "approved" or paste the observed shape summary</resume-signal>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| client (browser) -> Datagrok server | Spike issues real CRUD against the live server; throwaway entities created |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-04 | Denial of service | Spike test (project deletion on rollback) | mitigate | Spike wraps `projects.delete` in try/catch and logs error so operator can manually clean up if cascade fails — never silent partial-delete |
| T-15-SP-01 | Information disclosure | Spike test (probe outputs to log) | accept | Output may contain user friendly name / group names; this is by-design probe output, only the running operator sees it |
| T-15-SP-02 | Tampering | Spike test (creates throwaway umbrella Space) | mitigate | Umbrella Space name includes timestamp; cleanup via `spaces.delete` cascades to child Space and children (RESEARCH §"Nested Spaces Verification" point 4) |
</threat_model>

<verification>
- Spike test runs against live server, prints output covering all 8 enumeration points
- Operator confirms via checkpoint that the trim contract for Wave 1 is locked
- All throwaway entities cleaned up at end of spike (test self-cleans, but operator confirms via UI / `grok s projects list` if concerned)
- `src/package-test.ts` updated with the spike import line
</verification>

<success_criteria>
- `src/tests/publish-spike.ts` exists and runs without throwing
- `grok test --category Publishing-Spike` outputs enumeration for all 8 assumption probes
- Operator captures resolved shape and is unblocked to proceed to Wave 1
- No leftover entities in the live server post-spike (or operator has manually cleaned up via UI)
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md` when done with:
- Resolved values for Assumptions A1, A2, A3, A4, A6, A8 (one line each)
- Any unexpected shapes that affect Wave 1+ planning
- Decision: which Wave 1 plans (01, 02, 04) need locked-shape adjustments
</output>
