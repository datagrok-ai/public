---
phase: 15-read-only-publishing-foundation
plan: 03
type: execute
wave: 2
depends_on: ["15-00", "15-01"]
files_modified:
  - src/publishing/assert-published-shape.ts
autonomous: true
requirements: [PUB-03, PUB-06]
must_haves:
  truths:
    - "`assertPublishedShape(project)` reopens the saved Project from id, asserts the contract, throws on first violation with a specific error message"
    - "The contract checked: required `proteomics.published*` tags present (OR readable via belt-and-braces column), 7-column allowlist intact, PROTEIN_ID semType present on protein-id column, `df.name` matches publish-time name, volcano viewer present in `tv.viewers`, formula lines for FC threshold + p threshold present on the volcano viewer (PUB-06 / SC-2)"
    - "Helper is consumed by BOTH the test (Plan 08) AND defensively by the orchestrator (Plan 04 step 8) — single source of truth for the round-trip contract"
    - "On any contract violation, helper throws Error with a descriptive message including which assertion failed and what value was observed; volcano formula-line assertion message names a specific recovery hook so Plan 04's self-healing path (W-7) can detect and recover"
  artifacts:
    - path: "src/publishing/assert-published-shape.ts"
      provides: "Round-trip contract assertion"
      exports: ["assertPublishedShape", "PublishedShapeContract", "FORMULA_LINE_ASSERTION_PREFIX"]
  key_links:
    - from: "src/publishing/assert-published-shape.ts"
      to: "grok.dapi.projects.find"
      via: "reopen by id + .open()"
      pattern: "projects\\.find.*\\.open"
    - from: "src/publishing/assert-published-shape.ts"
      to: "src/publishing/publish-state.ts"
      via: "PUBLISHED_TAGS + META_COLUMNS + getPublishedMetadata imports"
      pattern: "import.*publish-state"
---

<objective>
Build the load-bearing round-trip gate for the entire phase. `assertPublishedShape(project, contract)` reopens the saved Project in a fresh session, asserts every observable invariant of the publish contract, throws with a specific error on first violation.

Purpose: Pitfall 3 (serializer-strip) is real (Phase 13 commit `e527d07ba1`). "Save succeeded" does NOT mean "reopen will look the same." The orchestrator (Plan 04 step 8) calls this helper INSIDE `publishAnalysis` before returning success — if reopen fails the contract, the Project is rolled back (or self-heals per W-7 / B-2). The test (Plan 08) calls the same helper to gate the regression suite.

**Revision B-2 addition:** PUB-06 / SC-2 mandates formula lines for FC and p thresholds. The Phase 13 evidence (commit `e527d07ba1`) is that look / filter config gets stripped by the serializer. Assertion 8a now explicitly checks for formula-line presence on the reopened volcano — and emits a distinctive prefix (`FORMULA_LINE_MISSING:`) so Plan 04's self-healing path can detect this specific failure mode and invoke the post-open recovery hook (B-2 + W-7) before rolling back.

Output: `src/publishing/assert-published-shape.ts` exporting `assertPublishedShape(project, contract)`. Single source of truth for what "successfully published" means.
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
@packages/Proteomics/CLAUDE.md
@packages/Proteomics/src/publishing/publish-state.ts
@packages/Proteomics/src/utils/proteomics-types.ts
@packages/Proteomics/src/viewers/volcano.ts
@packages/Bio/src/tests/projects-tests.ts

<interfaces>
<!-- Pulled from packages/Bio/src/tests/projects-tests.ts:49-73 (closest analog) -->

From @datagrok-libraries/test/src/test:
  awaitCheck(predicate: () => boolean, errorMessage: string, timeoutMs: number): Promise<void>
  delay(ms: number): Promise<void>

From datagrok-api:
  grok.dapi.projects.find(id: string): Promise<DG.Project>
  project.open(): Promise<void>
  grok.shell.tv: DG.TableView | null
  grok.shell.tv.dataFrame: DG.DataFrame
  grok.shell.tv.viewers: Iterable<DG.Viewer>
  v.type: string                     // DG.VIEWER.SCATTER_PLOT etc.
  v.getOptions(): {look?: {formulaLines?: string}; ...}  // formula lines live in look.formulaLines (JSON string of [{title, formula, color}])
  DG.VIEWER.SCATTER_PLOT = 'Scatter plot'

From Plan 01 (publish-state.ts):
  PUBLISHED_TAGS, META_COLUMNS, getPublishedMetadata(df), PublishedMetadata

From src/utils/proteomics-types.ts:
  SEMTYPE.PROTEIN_ID = 'Proteomics-ProteinId'
</interfaces>
</context>

<tasks>

<task type="auto" tdd="false">
  <name>Task 1: Define PublishedShapeContract + assertPublishedShape skeleton + FORMULA_LINE_ASSERTION_PREFIX constant</name>
  <files>src/publishing/assert-published-shape.ts</files>
  <read_first>
    - @packages/Bio/src/tests/projects-tests.ts (lines 49-73 — `dataFrameContainsColumns` + `checkViewerAdded` patterns — the only analog)
    - @packages/Proteomics/src/publishing/publish-state.ts (PUBLISHED_TAGS, META_COLUMNS, getPublishedMetadata, PublishedMetadata)
    - @.planning/phases/15-read-only-publishing-foundation/15-PATTERNS.md (Section 5 — analog and notes)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (success criterion 3 — exact contract)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"System Architecture Diagram" step 8 — non-negotiable gate #2)
    - @.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md (spike output — confirms which tags survive, drives reader path)
  </read_first>
  <action>
Create `src/publishing/assert-published-shape.ts`. Standard imports (`* as grok`, `* as DG`, `awaitCheck`, `delay` from `@datagrok-libraries/test/src/test`, `PUBLISHED_TAGS`, `META_COLUMNS`, `getPublishedMetadata`, `PublishedMetadata` from `./publish-state`, `SEMTYPE` from `../utils/proteomics-types`, `findColumn` from `../utils/column-detection`).

Export `FORMULA_LINE_ASSERTION_PREFIX = 'FORMULA_LINE_MISSING:' as const;` — this prefix is the load-bearing distinguishing marker that Plan 04 step 8's catch block uses to identify the formula-line-stripped failure mode and trigger the post-open recovery hook (B-2 / W-7).

Define `PublishedShapeContract` interface — captures what the publisher claimed at publish time, used as the expected baseline:
  - `expectedName: string` (the trimmed DF's name — `<source>_published_<date>`)
  - `expectedProjectName: string` (`Proteomics-Review-<slug>-v<N>-<date>`)
  - `expectedAllowlist: string[]` (the 7 column names the trim chose — proteinId, geneSymbol?, log2FC, pValue, adjPValue, significant, direction?)
  - `expectedMeta: PublishedMetadata`
  - `expectVolcano: boolean` (always true unless tested in isolation)
  - `expectEnrichment: boolean` (true iff D-05 carry happened)
  - `expectFormulaLines: boolean` (PUB-06 / SC-2 — always true when expectVolcano is true; settable to false only in unit-test isolation)

Define `assertPublishedShape(project: DG.Project, contract: PublishedShapeContract): Promise<void>` — async function, throws on first violation.

Skeleton:
- capture `projId = project.id`
- `grok.shell.closeAll(); await delay(100);`
- `const reopened = await grok.dapi.projects.find(projId); await reopened.open();`
- `await awaitCheck(() => !!grok.shell.tv?.dataFrame, 'assertPublishedShape: TableView never materialized after project.open()', 5000);`
- `const reDf = grok.shell.tv!.dataFrame;`
- Then run a series of `if (!cond) throw new Error('assertPublishedShape: <specific>: got <observed>, expected <expected>')` checks (Tasks 2, 3).
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/assert-published-shape.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "export interface PublishedShapeContract" src/publishing/assert-published-shape.ts | grep -q '^1$' &amp;&amp; grep -c "export async function assertPublishedShape" src/publishing/assert-published-shape.ts | grep -q '^1$' &amp;&amp; grep -c "FORMULA_LINE_ASSERTION_PREFIX" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -q "grok\.dapi\.projects\.find" src/publishing/assert-published-shape.ts &amp;&amp; grep -q "await reopened\.open" src/publishing/assert-published-shape.ts</automated>
  </verify>
  <done>File exists; `PublishedShapeContract` interface defined (now includes `expectFormulaLines`); `FORMULA_LINE_ASSERTION_PREFIX` constant exported; `assertPublishedShape` skeleton reopens by id, calls `await project.open()`, waits for shell.tv to materialize; strict TypeScript compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 2: Implement core assertions — name, tags+metadata column, semType, columns</name>
  <files>src/publishing/assert-published-shape.ts</files>
  <read_first>
    - @packages/Bio/src/tests/projects-tests.ts (lines 49-73 — column-presence assertion shape via `awaitCheck`)
    - @packages/Proteomics/src/publishing/publish-state.ts (META_COLUMNS keys — column-first read order)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md ("Belt-and-braces is the design philosophy" — column FIRST, tag SECOND)
    - @.planning/phases/15-read-only-publishing-foundation/15-00-SUMMARY.md (spike output — which tags survived clean save/reopen)
  </read_first>
  <action>
Add assertions to `assertPublishedShape` (after the reopen + DF-wait code in Task 1):

**ASSERTION 1 — df.name matches:** `if (reDf.name !== contract.expectedName) throw new Error('assertPublishedShape: df.name: got "' + reDf.name + '", expected "' + contract.expectedName + '"');`

**ASSERTION 2 — Project name matches (PUB-09):** compare `reopened.name !== contract.expectedProjectName`.

**ASSERTION 3 — Every allowlist column present:** loop `for (const colName of contract.expectedAllowlist)`. If `!reDf.col(colName)` throw with specific column name.

**ASSERTION 4 — Column count matches allowlist (no extras except hidden `_meta_*`):** filter out columns starting with `_meta_`; remaining count must equal `contract.expectedAllowlist.length`. On mismatch, throw an error listing the unexpected extra column names.

**ASSERTION 5 — PROTEIN_ID semType present on protein-id column:**
- Find col via `findColumn(reDf, SEMTYPE.PROTEIN_ID, ['protein id', 'majority protein id', 'accession'])`.
- If null or `col.semType !== SEMTYPE.PROTEIN_ID`, throw with hint: "Detectors may not have re-fired on reopen — re-assign in trim-dataframe.ts Step E."

**ASSERTION 6 — Required `proteomics.published*` metadata readable (tag OR column — Pitfall 3 mitigation):**
- Call `const meta = getPublishedMetadata(reDf);` (Plan 01 helper reads column FIRST, tag SECOND).
- If null: throw "getPublishedMetadata returned null — published flag not set or all critical fields unreadable".
- Then compare each required field against `contract.expectedMeta`: `target`, `deMethod`, `fcThreshold` and `pThreshold` (float epsilon `Math.abs(a-b) > 1e-9`), `version`, `publishId`, `includesEnrichment`, `publishedBy`. Each mismatch throws with field name + observed + expected.

**ASSERTION 7 — Belt-and-braces metadata column present (PUB-11):**
- For each key in `META_COLUMNS`: `if (!reDf.col(META_COLUMNS[k])) throw new Error('assertPublishedShape: metadata column missing: ' + META_COLUMNS[k]);`
- Catches the case where tag-only encoding silently passed (because tags survived in this cell) but the column was missing — PUB-11 specifically requires the belt-and-braces column.
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/assert-published-shape.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "getPublishedMetadata" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "META_COLUMNS\[" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; THROW_COUNT=$(grep -c "throw new Error" src/publishing/assert-published-shape.ts) &amp;&amp; [ "$THROW_COUNT" -ge 7 ]</automated>
  </verify>
  <done>At least seven explicit `throw new Error` assertion sites covering df.name, project name, each allowlist column, column count, PROTEIN_ID semType, getPublishedMetadata field equality, and every metadata column; descriptive messages naming the failed field; strict TypeScript compiles.</done>
</task>

<task type="auto" tdd="false">
  <name>Task 3: Implement volcano viewer + formula-line + enrichment cross-DF + supersede assertions</name>
  <files>src/publishing/assert-published-shape.ts</files>
  <read_first>
    - @packages/Bio/src/tests/projects-tests.ts (lines 65-73 — `checkViewerAdded` pattern via `awaitCheck`)
    - @packages/Proteomics/src/package.ts (volcano viewer add — search for `SCATTER_PLOT` usage)
    - @packages/Proteomics/src/viewers/volcano.ts (canonical volcano viewer factory — what props identify "this is the volcano"; how formula lines are added — typically `viewer.setOptions({look: {formulaLines: JSON.stringify([...])}})` or `viewer.props.formulaLines`)
    - @.planning/phases/15-read-only-publishing-foundation/15-RESEARCH.md (§"Pattern 1" — viewer present iteration; Phase 15-specific pitfall — volcano color binding survival)
    - @.planning/phases/15-read-only-publishing-foundation/15-CONTEXT.md (PUB-06 SC-2 — thresholds drawn as formula lines on reviewer first-paint)
  </read_first>
  <action>
Add to `assertPublishedShape` (after Assertion 7):

**ASSERTION 8 — Volcano viewer present (PUB-06):**
- Iterate `grok.shell.tv!.viewers` inside `awaitCheck` with 5000ms timeout.
- Predicate: any viewer with `v.type === DG.VIEWER.SCATTER_PLOT`. On found, capture into `volcanoViewer` outer variable for downstream use.
- Throw with hint "consider post-open recomputeVolcano in publish-project.ts" when `contract.expectVolcano` and not found.
- If spike output revealed that viewer props DID NOT survive reopen, skip prop check (defensive re-init fires in Plan 04). Otherwise, validate `getOptions()` carries the expected `x = log2FC`, `y = negLog10P` (or equivalent — match the `src/package.ts` precedent for "this is THE volcano").

**ASSERTION 8a — Volcano formula lines for FC + p thresholds (PUB-06 / SC-2 — REVISION B-2):**
- Only runs if `contract.expectFormulaLines === true` AND a volcano viewer was found in Assertion 8.
- Read formula lines from the volcano. Datagrok stores them on the viewer's look config. Use ONE of these access paths (executor picks the one that matches the installed Datagrok release; check `volcano.ts` for the precedent):
  - `volcanoViewer.getOptions().look?.formulaLines` (JSON string, parse it)
  - `(volcanoViewer.props as any).formulaLines` (also JSON string)
  - `(volcanoViewer as any).meta?.formulaLines`
- Parse the formula-line array. Each entry has shape `{title?, formula, color, ...}` where `formula` is a string like `'${y} = -log10(0.05)'` or `'${x} = 1'`.
- Required content (PUB-06):
  1. At least one formula line whose `formula` references the FC threshold value (`contract.expectedMeta.fcThreshold`). Match by substring: stringified `fcThreshold` appears in the formula expression (e.g. formula contains `'1'` when fcThreshold=1).
  2. At least one formula line whose `formula` references the p threshold value. Match by substring: stringified `pThreshold` (e.g. `'0.05'`) OR `-log10(pThreshold)` form.
- **THROW with `FORMULA_LINE_ASSERTION_PREFIX` prefix** — this is load-bearing for the self-healing path in Plan 04 step 8 catch block (W-7):
  ```
  throw new Error(FORMULA_LINE_ASSERTION_PREFIX + ' volcano reopened without formula lines for FC=' + contract.expectedMeta.fcThreshold + ' and/or p=' + contract.expectedMeta.pThreshold + '. The look/filter config was stripped by the serializer (Phase 13 e527d07ba1 evidence). Plan 04 step 8 catch should invoke the post-open recovery hook from tags proteomics.published_fc_threshold / proteomics.published_p_threshold and re-run this assertion once.');
  ```
- Skip this assertion entirely if `contract.expectFormulaLines === false` (allows unit-test isolation of the rest of the contract without the formula-line requirement).

**ASSERTION 9 — Enrichment DF presence (PUB-12, conditional):**
- If `contract.expectEnrichment === true`: `grok.shell.tables` must contain ≥ 2 DataFrames. Find the protein DF via its `_meta_published_includes_enrichment === 'true'` column, find the enrichment DF via its `proteomics.enrichment === 'true'` tag. If either missing, throw.
- Cross-DF wiring check is deferred to Plan 04 orchestrator + post-open recovery hook (W-5); this assertion only verifies presence.
- If `expectEnrichment === false`: skip.

**ASSERTION 10 — Project options bidirectional supersede pointers (PUB-10, conditional):**
- Per spike output A4 — if `Project.options` survival was verified, and `contract.expectedMeta.supersedes` is non-null:
  - Read `reopened.options['proteomics.supersedes']` (key per Plan 01 SUMMARY decision); throw on mismatch.
  - Same for `supersededBy` if `contract.expectedMeta.supersededBy` is non-null.
- Skip whichever side is null in the contract.
- Note for W-8: Plan 04 step 9 writes the `superseded_by` pointer to BOTH the prior Project's `options` AND the prior Project's DataFrame tag (belt-and-braces). This assertion checks the new project's `supersedes` — the prior's `superseded_by` is verified at the panel render time (Plan 06) and in Plan 08 Test 4.

End of function — no explicit return. Add JSDoc note: "Caller is responsible for `grok.shell.closeAll()` after the assertion if cleanup is desired. Plan 04 step 8 catches the FORMULA_LINE_ASSERTION_PREFIX-prefixed error specifically and invokes the post-open recovery hook before rolling back (W-7 self-healing path)."
  </action>
  <verify>
    <automated>cd packages/Proteomics &amp;&amp; npx tsc --noEmit src/publishing/assert-published-shape.ts 2>&amp;1 | grep -v '^$' | { ! grep -E 'error TS'; } &amp;&amp; grep -c "SCATTER_PLOT" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "expectEnrichment" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "expectVolcano" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "expectFormulaLines" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "FORMULA_LINE_ASSERTION_PREFIX" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "supersedes" src/publishing/assert-published-shape.ts | grep -qv '^0$' &amp;&amp; grep -c "formulaLines" src/publishing/assert-published-shape.ts | grep -qv '^0$'</automated>
  </verify>
  <done>Four additional assertions (8 volcano viewer presence with optional prop check, 8a volcano formula-line presence for FC + p thresholds throwing with FORMULA_LINE_ASSERTION_PREFIX-prefixed message, 9 optional enrichment-DF presence, 10 optional supersede pointer match); all use `awaitCheck` for async-materialization waits; strict TypeScript compiles.</done>
</task>

</tasks>

<threat_model>
## Trust Boundaries

| Boundary | Description |
|----------|-------------|
| publisher-claimed contract -> reopened actual shape | The contract is what we said we did; the reopened DF is what survived. Any mismatch = silent data loss in the published artifact |

## STRIDE Threat Register

| Threat ID | Category | Component | Disposition | Mitigation Plan |
|-----------|----------|-----------|-------------|-----------------|
| T-15-02 | Information disclosure (stale snapshot — silent partial publish) | assertPublishedShape | mitigate | Non-negotiable gate per CONTEXT.md "verify-and-rollback ... is one of two non-negotiable gates"; called inside publishAnalysis (Plan 04 step 8); failure rolls back via `dapi.projects.delete(project)` (or self-heals on FORMULA_LINE_MISSING per W-7) |
| T-15-SP-04 | Tampering (false positive — assertion passes but artifact wrong) | assertPublishedShape | mitigate | Plan 08 round-trip test deliberately fuzzes contract values to ensure assertions fire on violation; tested, not just trusted |
</threat_model>

<verification>
- TypeScript strict-mode compiles
- All 11 assertions present and reachable (every `throw new Error` describes the specific failed assertion; Assertion 8a uses FORMULA_LINE_ASSERTION_PREFIX)
- Function reopens via `dapi.projects.find(projId).open()` exactly once
- Uses `awaitCheck` with 5000ms timeout for async-materialization waits (Bio convention)
- `FORMULA_LINE_ASSERTION_PREFIX` is exported as a `const` and used inside the Assertion 8a throw
</verification>

<success_criteria>
- `src/publishing/assert-published-shape.ts` exports `assertPublishedShape` + `PublishedShapeContract` + `FORMULA_LINE_ASSERTION_PREFIX`
- Function is callable from both Plan 04 orchestrator and Plan 08 test (single-source-of-truth contract)
- Every assertion throws with field-name + observed + expected in the message
- Function does NOT mutate the source DF or close views beyond `grok.shell.closeAll()` at the start of reopen
- Assertion 8a throws with the FORMULA_LINE_ASSERTION_PREFIX so Plan 04 step 8 can detect this specific failure mode
</success_criteria>

<output>
Create `.planning/phases/15-read-only-publishing-foundation/15-03-SUMMARY.md` when done with:
- List of 11 assertions (now including 8a for formula lines)
- Notes on what spike output drove (e.g., "spike confirmed semType strip — Assertion 5 added re-assign hint message")
- Confirmation that this is the single source of truth (Plan 04 and Plan 08 will both call it)
- Confirmation that FORMULA_LINE_ASSERTION_PREFIX is exported for Plan 04 step 8 self-healing detection (B-2 / W-7)
</output>
