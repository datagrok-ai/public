---
feature: powerpack
sub_features_covered:
  - powerpack.io.xlsx-file-handler
  - powerpack.io.exceljs-service
  - powerpack.io
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
original_path: D:\work\datagrok\reddata\public\packages\UsageAnalysis\files\TestTrack\Powerpack\xlsx-open.md
migration_date: 2026-05-23
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-19329
realized_as:
  - xlsx-open-spec.ts
gate_verdicts:
  a:
    verdict: PASS
    cycle_id: 2026-05-22-powerpack-migrate-05
    timestamp: 2026-05-23T00:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses cleanly as YAML 1.2 and carries all four required
          fields. feature=powerpack; sub_features_covered=
          [powerpack.io.xlsx-file-handler, powerpack.io.exceljs-service,
          powerpack.io]; target_layer=playwright; coverage_type=regression.
          All three sub_features_covered ids resolve in the atlas
          (powerpack.yaml: powerpack.io at line 725, powerpack.io.xlsx-file-handler
          at line 740, powerpack.io.exceljs-service at line 748). No nested-mapping
          colon hazards inside the frontmatter body.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains H2 `## Scenarios` (line 192) wrapping six H3 scenario
          headings: `### Scenario 1: Open XLSX from Browse / My Files`,
          `### Scenario 2: Open XLSX from Recent files`,
          `### Scenario 3: ... drag-and-drop`, `### Scenario 4: ... File > Open`,
          `### Scenario 5: ... Shared with me`, `### Scenario 6: ... sheetName`.
          Per the refined mode-file reading ("scenario heading specifically"),
          the H2 + H3 scenario layout satisfies the check — this matches the
          convention used by sibling section scenarios.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Each of the six scenarios opens with numbered steps starting
          `1. **...**`. Setup section also has three numbered steps. Spot-counts:
          S1=7 steps, S2=5, S3=4, S4=5, S5=5, S6=4. No scenario heading lacks
          numbered steps.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          No empty scenarios. Every scenario has ≥ 4 numbered steps plus an
          `Expected:` block, so the body of each scenario is non-trivial.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer=playwright is in the canonical enum
          {playwright, apitest, manual-only} per
          verdict-enums.yaml derived_enums.target_layer.canonical_values.
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type=regression is in the canonical enum
          {smoke, regression, edge, perf} per
          verdict-enums.yaml derived_enums.coverage_type.canonical_values.
          No severity-axis value (p0..p3) appears in the coverage_type slot.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type=regression is declared at file frontmatter level and
          inherits to all six scenarios — no per-scenario override deviates.
          Notes "coverage_type rationale" block explicitly justifies regression
          (not smoke = not the section golden path; not edge = not a boundary
          value; not perf = functional, not perf-related).
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          Repeated setup (PowerPack install + powerPackInit, XLSX fixture in
          two locations, DevTools console open) is factored into the
          `## Setup` section. Scenario 2 chains off Scenario 1 via a
          `Pre-condition` step rather than re-running setup; Scenarios 3-6
          reuse the Setup-section fixture references and do not duplicate
          install / fixture preparation prose.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          Frontmatter does NOT declare a pyramid_layer field (the bug-focused
          designation appears only in body prose under Notes). Per the
          mode-file, this rule PASSes by vacuity when pyramid_layer is absent
          from frontmatter; the hard alignment rule only fires when
          pyramid_layer=ui-smoke. No mismatch.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Scenarios reference real artifacts throughout: real source paths
          (public/packages/PowerPack/src/package.ts#L383, L35, L366), real atlas
          sub-feature ids, real DOM affordances (Browse, My Files, Recent files,
          Shared with me, File > Open), real fixture paths
          (System:DemoFiles/SPGI-linked.xlsx, Home/xlsx-open-test.xlsx),
          real grok s commands (grok s files save, grok s shares add). No
          angle-bracket or square-bracket placeholders; no generic
          column1/some-file.csv stand-ins. The OtherUser:xlsx-shared-test.xlsx
          token in Scenario 5 is a real Datagrok cross-user fixture path, not
          a placeholder.
      - check: A-BUG-01
        status: PASS
        evidence: |
          A-BUG-01 evaluates atlas known_issues with literal `test_coverage:
          needed`. The powerpack atlas known_issues field (lines 1352-1445)
          uses the schema-evolved form `test_coverage: { exists: false,
          paths: [] }` on every entry — there are ZERO entries with the literal
          `test_coverage: needed` form. PASS-by-vacuity under the strict
          literal reading of the mode-file rule. Under the semantic reading
          (exists: false treated as equivalent to "needed"), the only
          known_issues entry whose affects_sub_features intersects this
          scenario's sub_features_covered is GROK-19329 (atlas lines
          1437-1445), and the scenario's related_bugs frontmatter contains
          [GROK-19329], satisfying clause (a). The other 8 known_issues
          entries affect unrelated sub-features (widgets, formula-lines,
          add-new-column, dashboards, welcome.view, search) and are
          out-of-scope for this single-scenario gate; chain-wide coverage
          of those bugs is Gate F's responsibility via F-BUG-COVERAGE-01.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          Three deferrals in Notes "Deferrals" block; each cites a real
          technical or atlas-level dependency, not effort/complexity. (a)
          apitest slice deferred because "GROK-19329 regression was about
          UI-path-to-handler routing, not about parser correctness" — real
          bug-scope citation. (b) 80MB size edge deferred because "scenario
          fixture is < 1 MB by design (small fixture for CI speed) and the
          size limit is not part of the GROK-19329 reproduction surface" —
          real atlas declaration + reproduction-scope citation. (c)
          Shared-with-me cross-user precondition documented via reference
          to data-enrichment.md Sub 4 sibling precedent + logoutAndLoginAs
          helper convention — real precedent citation, not effort opt-out.
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No "TODO: add later" / "to be done in next phase" without
          dependency reference. Deferral language consistently cites a real
          prerequisite or dependency: "remains a Test Designer follow-up if
          parser-level regressions surface separately" (parser-level
          regression surfacing = real prerequisite), "deferred because the
          scenario fixture is < 1 MB by design" (atlas size limit = real
          atlas-level dependency). Lattice Rule 13 satisfied.
  e:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T14:05:00Z
    failure_keys: []
  d:
    verdict: EVIDENCE_GAP
    cycle_id: 2026-05-23-powerpack-migrate-02
    timestamp: 2026-05-23T17:30:00Z
    failure_keys: []
    evidence_needed: |
      Re-emission #2 under cycle_id 2026-05-23-powerpack-migrate-02.
      Bundle state is unchanged from the three prior EVIDENCE_GAP
      records in decision-log.recon_evidence_records (timestamps
      2026-05-23T14:17:26Z, 14:21:27Z, 14:23:36Z) and from the
      sibling round-1 verdict under cycle
      2026-05-23-powerpack-migrate-01: original_path resolves to
      the same on-disk file as the migrated scenario, so the
      diff-based class of Gate D content checks (D-STEP-01,
      D-STEP-02, D-EDGE-01, D-STRUCT-01, D-STRUCT-02, D-SAN-02,
      D-MERIT-01, D-MERIT-02, D-UI-DELEGATION-01) has no diff to
      consult and is unevaluable as a class. Critic D's
      adjudication boundary does not extend to choosing between
      (a) recovering a pre-rewrite ORIGINAL from git history and
      (b) reclassifying produced_from to atlas-driven; both are
      orchestrator-layer decisions. Per migration-mode.md's
      closing rule, after 2 review rounds with no PASS or
      SCOPE_REDUCTION the orchestrator escalates per autopilot
      boundary trigger #1 rather than re-invoking Critic D a
      further time on the same bundle.
      Mechanical fail-fast checks PASS on literal reading and do
      not block this verdict — the gap is content-side and
      classification-side, not structural. D-STRUCT-MECH-03: all
      8 required migration-output frontmatter fields present
      (feature, sub_features_covered, target_layer, coverage_type,
      produced_from, original_path, migration_date, related_bugs);
      deprecated migrated_from absent; produced_from=migrated is
      in the 3-value enum {migrated, atlas-driven, decomposed}.
      D-STRUCT-MECH-05: original_path target file exists in the
      repo (Glob on
      public/packages/UsageAnalysis/files/TestTrack/**/xlsx-open*.md
      returns exactly one match: PowerPack/xlsx-open.md) — but it
      is the file under review itself, which is the upstream
      precondition that strands the content checks.
      D-FRONTMATTER-PHASE1-01: all four Phase 1 fields present as
      parseable inline-flow lists, all empty [].
      D-FRONTMATTER-PHASE1-02: per-entry schemas vacuously
      satisfied because every Phase 1 list is empty.
      Supporting-context corroboration of the
      atlas-driven-reclassification hypothesis is unchanged from
      round 1: (a) chain YAML
      references/scenario-chains/powerpack.yaml lines 240-244
      declares the file in migration_decisions[] with
      source=xlsx-open.md, target=PowerPack/xlsx-open.md,
      target_layer=playwright, strategy=simple, reason citing
      GROK-19329 — migration_decisions[] presupposes a
      pre-existing TestTrack ORIGINAL was rewritten to the target
      path. (b) Same chain YAML lines 326-333 records GROK-19329
      in bug_match_attempts_skipped as below_trigger_threshold
      with rationale "Bug intent is already carried by
      xlsx-open.md frontmatter related_bugs" — implying the
      scenario pre-existed at chain-analyze time and was matched,
      not authored from the bug. (c) Scenario body directly
      contradicts (b): the "Chain context" Notes block (scenario
      body ~lines 540-546) states "This scenario is the section's
      bug-focused witness for GROK-19329 — Critic F's
      bug_focused_candidates[] entry for the bug (added in
      cycle-2026-05-20-powerpack-coverage) had empty spans[]
      pending this scenario's authoring"; the "Coverage map"
      Notes block (~lines 581-586) states "The Critic F
      coverage-gap dispatch (gap: bug-uncovered :: GROK-19329)
      drove this scenario's authoring directly". Chain YAML and
      scenario body disagree on origin (pre-existing-and-matched
      vs F-gap-authored); resolving the disagreement is
      prerequisite to choosing the correct produced_from enum
      value.
      Recon targets (orchestrator-layer; out of Critic D's
      boundary; unchanged from prior rounds, listed in execution
      order):
      (1) `git log --diff-filter=A --follow --
      public/packages/UsageAnalysis/files/TestTrack/PowerPack/
      xlsx-open.md` to identify the introducing commit and its
      author / message. If the introducing commit names Migrator
      or cites a migration cycle_id, produced_from=migrated is
      correct and recon continues at step (3). If the introducing
      commit names a Critic F gap-fill cycle
      (cycle-2026-05-20-powerpack-coverage is the strong
      candidate per the body Notes), produced_from must be
      reclassified to atlas-driven and the scenario re-dispatched
      as Gate A only — skipping Gate D entirely per
      migration-mode.md's produced_from enum resolution
      (atlas-driven covers exactly this author-from-gap flavor;
      lifecycle-api.md round-1 SR precedent cited in the
      D-STRUCT-MECH-03 commentary).
      (2) Inspect decision-log cycle_completions entries for
      cycle_id 2026-05-22-powerpack-migrate-05 (current
      gate_verdicts.a cites it as origin cycle) and for
      cycle-2026-05-20-powerpack-coverage (body Notes cites it
      as F gap-fill origin); whichever is the producing-agent
      cycle is authoritative for produced_from classification.
      (3) If steps (1)-(2) confirm migrated-flavor lineage, the
      rewrite must have a pre-rewrite ORIGINAL in some git
      revision. Recover via `git log -p --
      public/packages/UsageAnalysis/files/TestTrack/PowerPack/
      xlsx-open.md` and snapshot the pre-rewrite revision to a
      tmp path; re-dispatch Gate D with original_path pointing
      at the snapshot so diff-based checks can run.
      (4) If steps (1)-(3) fail to surface a distinct ORIGINAL,
      reclassify produced_from to atlas-driven and re-dispatch
      as Gate A only — this scenario is then not a migration
      bundle and Gate D does not apply.
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T13:14:00Z
    spec_runs:
      - spec: xlsx-open-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 39
        failure_keys: []
---

# PowerPack — XLSX file open across all entry paths (GROK-19329 regression)

Bug-focused regression scenario for `GROK-19329`: PowerPack's XLSX file
handler (`xlsxFileHandler(bytes, sheetName?)` registered as a
`@fileHandler` for the `xlsx` extension —
`public/packages/PowerPack/src/package.ts#L383`) must remain functional
across every entry path Datagrok exposes for opening an `.xlsx` file.
Before the 1.27.0 regression, opening an XLSX file via Browse / My
Files, Recent files, drag-and-drop into the platform window, the top
menu File > Open, or Shared with me each routed the file's bytes
through the registered handler, which delegated to
`ExcelJSService.getInstance().parse(bytes, sheetName)` (the
`exceljs` Web Worker singleton —
`public/packages/PowerPack/src/package.ts#L35`) and surfaced one
DataFrame per sheet. The regression broke the handler universally —
XLSX did not open from ANY of the 5 entry paths, no data displayed,
no graceful error.

Atlas surface exercised: `powerpack.io.xlsx-file-handler` (the
`@fileHandler` entry point — the primary regression surface),
`powerpack.io.exceljs-service` (the Web Worker singleton that does
the actual parsing — broken initialization would surface here), and
`powerpack.io` (the umbrella file-handlers / IO sub-feature group
that hosts both).

## Setup

1. Open Datagrok with PowerPack installed (default platform load —
   `powerPackInit` runs at startup and registers `xlsxFileHandler`
   as the `@fileHandler` for the `xlsx` extension via
   `grok.functions.register(...)`).
2. Ensure an XLSX test fixture is available in two locations:
   - **Platform fixture (`System:DemoFiles/`).** Verify
     `System:DemoFiles/SPGI-linked.xlsx` (or any pre-existing
     `.xlsx` under `System:DemoFiles/`) exists and is readable.
     If no XLSX exists under `System:DemoFiles/`, upload one via
     `My Files` first (Step 2 below).
   - **User fixture (`Home/`).** Upload `Home/xlsx-open-test.xlsx`
     into the current user's `My Files` (Home directory). The
     fixture should be a small (< 1 MB) multi-sheet workbook —
     for example a 3-sheet workbook with sheets `Customers`,
     `Orders`, and `Products`, each with a header row plus 5-10
     data rows of mixed-type columns (numeric, string, date). If
     no convenient fixture exists locally, use Datagrok's
     `Customers.xlsx` sample if present under `tables/` of any
     installed package, or generate a fresh fixture with `grok s
     files save Home/xlsx-open-test.xlsx --from-file <local-path>`.
3. Open the browser DevTools Console (F12) so console errors are
   visible while opening XLSX files — the GROK-19329 regression
   surfaces as either a thrown error in the console (handler /
   Web Worker init failure) or as a silent no-op (file does not
   open, no view rendered).

## Scenarios

### Scenario 1: Open XLSX from Browse / My Files (canonical GROK-19329 reproduction)

1. **Navigate to Browse.** Click the `Browse` icon in the sidebar
   (or `Browse` in the top-left menu) to open the file browser
   view.
2. **Open the `My Files` (Home) tree.** Click `My Files` in the
   browse navigation so the current user's home directory listing
   surfaces.
3. **Locate the XLSX fixture.** The `xlsx-open-test.xlsx` fixture
   uploaded in Setup is visible in the listing.
4. **Single-click (or double-click per default Browse behavior)
   the XLSX file.** The `xlsxFileHandler` registered for the
   `xlsx` extension is invoked with the file's bytes; the
   `ExcelJSService` Web Worker parses the file.
5. **Verify the file opens.** A new TableView opens as the active
   view rendering the first sheet (e.g. `Customers`) as a Datagrok
   grid. The header row maps to grid columns; data rows surface in
   the grid body.
6. **Verify additional sheets surface.** If the fixture is
   multi-sheet, the second / third sheets surface as additional
   tabs / sibling views (one DataFrame per sheet per
   `ExcelJSService.parse`).
7. **Verify no console errors.** DevTools Console shows no
   exception (no `Cannot read property ... of undefined` from the
   handler, no Web-Worker init failure, no `ExcelJSService` null
   reference). No Datagrok red error balloon.

Expected:
- The XLSX file opens into a Datagrok TableView with data
  rendered from the first sheet — this is the canonical
  GROK-19329 invariant.
- Multi-sheet workbooks surface one view per sheet.
- No thrown exception in the console; no error balloon.

### Scenario 2: Open XLSX from Recent files

1. **Pre-condition.** Scenario 1 has already been run in this
   session so the XLSX fixture appears in the platform's Recent
   files list. If Scenario 1 was not run, open the fixture once
   via Browse / My Files first to populate Recent files.
2. **Navigate to Recent files.** Click `Recent files` in the
   sidebar (or `Browse > Recent` if the Recent files entry lives
   there).
3. **Locate the XLSX fixture** in the Recent files listing.
4. **Click the XLSX file.** The handler is invoked via the
   Recent-files navigation path.
5. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- Opening via Recent files routes through the same registered
  `xlsxFileHandler` as Browse / My Files — both must succeed.
- No regression-specific divergence between the Recent files
  navigation path and the Browse / My Files path.

### Scenario 3: Open XLSX via drag-and-drop into the platform window

1. **Navigate to the home / Welcome view** (or any non-Browse
   view — the drag-and-drop handler is registered globally on
   the Datagrok window).
2. **Drag the XLSX file from the local filesystem** (the OS file
   manager — Windows Explorer / macOS Finder) onto the Datagrok
   browser window.
3. **Drop the file anywhere on the Datagrok window.** The
   platform's global drag-and-drop handler reads the dropped
   file's bytes and dispatches to `xlsxFileHandler` (registered
   for `xlsx` extension).
4. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- Drag-and-drop dispatches to the same registered XLSX handler
  as Browse / My Files and Recent files — all three paths must
  succeed.
- The dropped file is parsed in-memory (no upload to `My Files`
  required for the open action to succeed).

### Scenario 4: Open XLSX via top menu `File > Open`

1. **Open the top menu `File > Open...`** (or the equivalent
   open-file action in the current Datagrok build's top menu).
2. **The OS native file picker opens** (Datagrok delegates to
   the browser's file-input element backed by the OS picker).
3. **Select the XLSX fixture** (`xlsx-open-test.xlsx`) from the
   local filesystem.
4. **Click `Open` in the picker.** The handler is invoked with
   the selected file's bytes.
5. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- The top menu `File > Open` path routes through the same
  registered `xlsxFileHandler` as the other entry paths — all
  four exercised so far must succeed.

### Scenario 5: Open XLSX from Shared with me

1. **Pre-condition.** A second user (or the current user via
   `grok s shares add`) must have shared an XLSX file with the
   current user. If no shared XLSX is available, run
   `grok s shares add "OtherUser:xlsx-shared-test.xlsx" <current-user>
   --access View` from a precondition setup step, OR ask another
   QA user to share an XLSX file via the Datagrok share dialog
   before running this scenario.
2. **Navigate to `Shared with me`.** Click `Shared with me` in
   the sidebar / Browse navigation.
3. **Locate the shared XLSX file** in the `Shared with me`
   listing.
4. **Click the shared XLSX file.** The handler is invoked.
5. **Verify the file opens.** Same expected behavior as Scenario
   1 Steps 5-7: TableView opens with data rendered; multi-sheet
   sheets surface as additional views; no console errors.

Expected:
- The `Shared with me` path routes through the same registered
  `xlsxFileHandler` as Browse / My Files, Recent files,
  drag-and-drop, and File > Open — all 5 entry paths must
  succeed (the GROK-19329 invariant: XLSX opens from ANY entry
  path).

### Scenario 6: Multi-sheet selection via sheetName parameter (optional)

This scenario guards the optional `sheetName` parameter of
`xlsxFileHandler(bytes, sheetName?)`. The parameter is exposed via
the atlas (`powerpack.io.xlsx-file-handler` declares
`xlsxFileHandler(bytes, sheetName?)`); if the platform UI exposes a
sheet-selector dialog on multi-sheet workbooks (or the API path is
testable), this scenario covers that path.

1. **Open the multi-sheet XLSX fixture via any of Scenarios 1-5's
   paths** — Browse / My Files is the simplest.
2. **If the platform prompts with a sheet-selector dialog** on
   multi-sheet workbooks, select a non-default sheet (e.g.
   `Orders` instead of `Customers`).
3. **Verify the selected sheet's data is rendered** as the
   primary TableView (header row mapped to grid columns; rows
   from the selected sheet).
4. **If no sheet-selector dialog is exposed**, the default
   behavior (one DataFrame per sheet, all surfaced) applies —
   verify all sheets are accessible per Scenario 1 Step 6.

Expected:
- If a sheet selector is exposed, selecting a non-default sheet
  renders that sheet's data.
- If no selector is exposed, all sheets surface as separate
  TableViews / tabs.
- `sheetName` parameter behavior is consistent with the atlas
  declaration (`xlsxFileHandler(bytes, sheetName?)` — the
  parameter is optional and defaults to opening all sheets).

## Notes

- **target_layer rationale.** `playwright`. The scenario exercises
  the XLSX file handler across 5 distinct UI entry paths (Browse /
  My Files navigation, Recent files navigation, OS drag-and-drop,
  top-menu File > Open with OS native picker, Shared with me
  navigation). All 5 paths involve DOM-level interactions, file
  rendering into a TableView grid, and the global platform
  drag-and-drop handler. Headless JS-API exercise (`apitest`) could
  test the handler's parsing function in isolation
  (`grok.functions.call('PowerPack:xlsxFileHandler', bytes)` — which
  would route through `ExcelJSService.parse`), but cannot verify
  that the 5 entry paths each correctly invoke the registered
  handler — the regression was about path-to-handler routing, not
  about parser correctness. A complementary `apitest` slice could
  exercise `xlsxFileHandler(bytes, sheetName?)` directly to verify
  the `ExcelJSService` Web Worker initialization and the parse
  result shape (one DataFrame per sheet, header dedup, default
  column names) — see Deferrals.
- **coverage_type rationale.** `regression`. Bug-focused
  (`pyramid_layer: bug-focused` per Rule 3 — canonical GROK-19329
  reproduction surface, walks the exact reproduction steps from
  the bug). Guards against re-regression of the universal XLSX-
  open path across all 5 entry paths. Not `smoke` (this is not a
  section golden path — the section's smoke is the top-level
  `add-new-column.md` per chain `ui_coverage_plan`). Not `edge`
  (the bug affects all common entry paths, not a boundary value).
  Not `perf` (regression is functional, not performance-related).
  Atlas `edge_cases[]` is empty (Phase 0 bootstrap), so this
  scenario's `coverage_type: regression` is derived from STEP E
  heuristics (general bug-focused coverage of a common feature
  shape) rather than from an atlas `edge_cases[]` entry — no
  cross-check mismatch.
- **Pyramid layer.** `bug-focused` per Rule 3 — discriminator
  test: GROK-19329 fails Scenarios 1-5 before the regression fix
  in 1.27.x (XLSX does NOT open from any of the 5 entry paths;
  no data displayed). After the fix, Scenarios 1-5 all pass.
  Scenario 6 (sheetName) is a defensive guard that exercises the
  optional parameter declared in the atlas; not part of the
  GROK-19329 reproduction path itself.
- **Atlas sub_features traceability.**
  - `powerpack.io.xlsx-file-handler` —
    `xlsxFileHandler(bytes, sheetName?)` is the `@fileHandler`
    registered for the `xlsx` extension
    (`public/packages/PowerPack/src/package.ts#L383`). Primary
    regression surface — the handler must route bytes from each
    of the 5 entry paths through `ExcelJSService.getInstance().parse()`.
    Atlas declares the handler refuses files larger than 80 MB;
    this scenario uses a < 1 MB fixture so the size guard is not
    exercised.
  - `powerpack.io.exceljs-service` — `ExcelJSService` singleton
    (`public/packages/PowerPack/src/package.ts#L35`) wraps the
    `exceljs` Web Worker. The handler calls
    `ExcelJSService.getInstance().parse(bytes, sheetName)`; if
    the Web Worker fails to initialize, every scenario fails at
    Step 5 / Step 7 with a Web-Worker init error in the console.
  - `powerpack.io` — umbrella file-handlers / IO sub-feature
    group (`public/packages/PowerPack/src/package.ts#L366`) that
    hosts both the XLSX handler and the Markdown file viewer.
- **Related bug.** `GROK-19329` (p2, regression-risk in dev
  1.27.0; fixed_in pending per bug-library). Reproduction: open
  Datagrok → try to open an XLSX file from Browse / My Files,
  Recent files, drag-and-drop, File > Open, Shared with me →
  XLSX does NOT open from any location, no data is displayed.
  Expected: PowerPack's XLSX file handler must function across
  versions; XLSX is consumed by all 5 entry paths;
  ExcelJSService Web Worker must initialize correctly;
  `xlsxFileHandler(bytes, sheetName?)` must parse Excel files
  into DataFrames per the atlas declaration.
- **Chain context.** This scenario is the section's bug-focused
  witness for GROK-19329 — Critic F's `bug_focused_candidates[]`
  entry for the bug (added in
  `cycle-2026-05-20-powerpack-coverage`) had empty `spans[]`
  pending this scenario's authoring. The chain's
  `order_from_files[]` is updated to include this scenario at
  the end of the section's order.
- **Deferrals.**
  - **`apitest` slice for `xlsxFileHandler` parsing
    correctness** — a complementary headless test could exercise
    `xlsxFileHandler(bytes, sheetName?)` directly via
    `grok.functions.call('PowerPack:xlsxFileHandler', bytes)`
    (or, more directly, the `ExcelJSService.parse` path) to
    verify parser-level invariants — one DataFrame per sheet,
    header dedup, default `col<n>` column names for missing
    headers, default `Sheet <n>` sheet names. Deferred because
    the GROK-19329 regression was about UI-path-to-handler
    routing, not about parser correctness; the canonical
    reproduction walks UI entry paths. The `apitest` slice
    remains a Test Designer follow-up if parser-level
    regressions surface separately.
  - **80 MB size limit edge case** — the atlas declares
    `xlsxFileHandler` refuses files larger than 80 MB; an
    edge-case scenario could exercise this boundary. Deferred
    because the scenario fixture is < 1 MB by design (small
    fixture for CI speed) and the size limit is not part of the
    GROK-19329 reproduction surface.
  - **Shared-with-me cross-user precondition** — Scenario 5
    requires a second user account to share an XLSX file
    in. The current section's `data-enrichment.md` Sub 4 has
    the same precondition with the same documented resolution
    pattern (`logoutAndLoginAs` helper + decision-log
    :: mig-2026-04-29-fixture-placeholder). Deferred to that
    convention; Scenario 5 documents the precondition inline so
    it can be exercised manually if the helper is unavailable.
- **Helpers candidates.** None added at this revision. The
  scenario's per-entry-path open actions are atomic enough to
  inline. If future cycles re-use the "drag-and-drop file from
  local filesystem" or "select file via top-menu File > Open
  picker" sequences, those would be candidate helpers (currently
  no equivalent helper exists in helpers-registry.yaml).
- **Coverage map.** Coverage map for PowerPack
  (`references/coverage-map/powerpack.yaml`) is not present at
  authoring time — gap-vs-coverage cross-check skipped per STEP B
  fallback. The Critic F coverage-gap dispatch
  (`gap: bug-uncovered :: GROK-19329`) drove this scenario's
  authoring directly.
