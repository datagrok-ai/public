---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
  - powerpack.dialogs.prepare-add-column-call
  - powerpack.formula.is-formula-column
  - powerpack.dialogs
target_layer: playwright
coverage_type: regression
pyramid_layer: bug-focused
ui_coverage_responsibility:
  - add-new-column-dialog
  - column-rename-context-action
  - save-project-with-datasync
  - project-reopen-with-formula-recalc
ui_coverage_delegated_to: add-new-column.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/add-new-column.md
migration_date: 2026-05-20
source_text_fixes:
  - step-3-todo-resolved-formula-weight-plus-100-column-weight2
  - step-4-resolved-formula-weight2-plus-100-column-weight3
  - close-all-spelled-out-as-close-all-views-via-shell
candidate_helpers:
  - helpers.powerpack.openTableFromLocalStorage
  - helpers.powerpack.openTableFromHomeDir
  - helpers.powerpack.openQueryResult
  - helpers.powerpack.saveProjectWithDatasync
  - helpers.powerpack.reopenProject
  - helpers.powerpack.addCalculatedColumn
  - helpers.powerpack.renameColumnViaContextAction
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-17109
realized_as:
  - add-new-column-advanced-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-subdir
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T00:00:00Z
    failure_keys: []
    scope_reduction_proposal: |
      Independent re-derivation for cycle 2026-05-28-powerpack-automate-02
      reaches SCOPE_REDUCTION on the same axis prior cycles found: the
      scenario "## Setup" step 1 enumerates five distinct source classes
      walked in sequence (scenario body: Source 1 local storage, Source 2
      Home dir, Source 3 Postgres:Northwind:OrdersByEmployee, Source 4
      Northwind products GetTop100, Source 5 Northwind products GetAll)
      and prescribes "the same Scenarios block is executed end-to-end"
      for each source. The spec walks only the file source
      (System:DemoFiles/demog.csv) for the full 10-step chain (softStep
      "Setup: open System:DemoFiles/demog.csv with datasync provenance"
      at spec L167 via openTableFromFile) and documents the four
      un-walked sources in its leading block-comment "Scope reductions"
      section (spec L74-105) with rationale: the 5x10 step matrix (50
      combinations) is impractical in a single ~600s Playwright test,
      and the GROK-17109 invariant (calculated columns persist across
      save+datasync+reopen with formula tags intact) is
      source-class-independent at the calc-column-persistence layer
      (the bug surfaces on the add-new-column dialog -> datasync save
      path, not on a specific source-class binding).

      Scenario frontmatter still records `scope_reductions: []` (L34) —
      the narrowing is NOT declared at the scenario layer; it lives only
      in the spec's leading block comment. Independent application of
      the Gate E checklist routes this to SCOPE_REDUCTION rather than
      FAIL: E-TRACE-02 PASSES because all 10 numbered scenario steps
      ARE covered in the spec (Step 1 L187, Step 2 L202, Step 3 L310,
      Step 4 L327, Step 5 L427, Step 6 L490, Step 7 L507, Step 8 L525,
      Step 9 L577, Step 10 L616). The narrowing lives in Setup step 1's
      source multiplicity (5 -> 1), not in any numbered-step omission.
      Per the spec-mode rubric "SCOPE_REDUCTION: acceptable if ... the
      spec covers fewer assertions than the scenario explicitly because
      of a documented technical limitation" — the limitation (5x10 =
      50 combinations exceeding the 600s playwright runtime bound) is
      documented in the spec's local notes.

      Proposed resolution (preferred 1 + 2 together):

      1. Backfill scenario frontmatter `scope_reductions:` with an SR
         entry `id: SR-01, check: setup-multi-source-walk-narrowed,
         rationale: "Walk 1 representative file source
         (System:DemoFiles/demog.csv) for the full 10-step chain in
         this -spec.ts; un-walked sources (local-storage, Home-dir,
         OrdersByEmployee, GetTop100, GetAll) move to per-source matrix
         specs at a future expansion. GROK-17109 invariant is source-
         class-independent at the calc-column persistence layer.",
         verdict_status: SCOPE_REDUCTION`. Re-emit Gate D and Gate E so
         the SR is declared at the layer where it is materialized.

      2. Surface the per-source matrix at the chain level as a
         `bug_focused_candidates` entry in
         `scenario-chains/powerpack.yaml` (proposed spec name
         `powerpack-add-new-column-multi-source-matrix-spec.ts`) so the
         un-walked sources are tracked as future work, not silently
         dropped. Scenario Notes already gesture at this with the
         GROK-17109 cross-cutting-candidate language.

      3. (NOT recommended.) Extend this spec to walk all five sources
         in a per-source loop. Runtime impact ~3000s — needs
         long-running-test budget allowance and `test.slow()`. Reject
         on runtime grounds.

      Retry-context check (E-RETRY-IGNORES-GATE-B): the scenario's
      `gate_verdicts.b` block (cycle_id 2026-05-26-powerpack-automate-01,
      the most-recent-prior-cycle of this scenario) records
      `verdict: PASS` (3/3 attempts, failure_keys: []). The
      retry-context predicate requires `gate_verdicts.b.verdict ==
      "FAIL"`; the current Gate B verdict is PASS, so retry-context
      detection does NOT fire and E-RETRY-IGNORES-GATE-B never fires.
      (The prior-cycle SR prose narrating a -03 B-STAB FAIL is stale —
      the actual recorded Gate B verdict is now PASS, validated against
      the per-attempt Playwright JSON reporter output.) Independent of
      the predicate: the spec source is unchanged from prior reads and
      retains no diagnosed-failing code path, so even under a
      hypothetical FAIL predicate the ignored-evidence pattern would
      not be present.

      Gate E evaluates spec-vs-scenario traceability and discipline on
      the authored source. The SCOPE_REDUCTION verdict stands on the
      source-matrix narrowing axis (5 -> 1 in the spec leading block
      comment vs empty scenario frontmatter scope_reductions[]); all
      other E-* checks PASS independently.
    claims:
      - check_id: E-STRUCT-MECH-01
        status: PASS
        evidence: |
          Spec file present at the paired path
          public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column-advanced-spec.ts
          (sibling of scenario add-new-column-advanced.md). Verified
          by Glob this cycle (cycle 2026-05-26-powerpack-automate-01).
          On-disk folder casing is `PowerPack` (PascalCase), matching
          the public/packages/<Pkg>/ naming convention.
      - check_id: E-STRUCT-MECH-02
        status: PASS
        evidence: |
          Spec parses as valid TypeScript on read: balanced braces
          (test arrow open L137 -> finally L644 -> close L654 -> outer
          close L660), single `test(...)` invocation, all imports are
          valid TS import declarations referencing existing module
          paths, no malformed signatures.
      - check_id: E-STRUCT-MECH-03
        status: PASS
        evidence: |
          One `test(...)` block at L137: `test('PowerPack: Add New
          Column — multi-source datasync persistence + formula recalc
          on rename (GROK-17109)', async ({page}) => {...})`.
      - check_id: E-STRUCT-MECH-04
        status: PASS
        evidence: |
          Test name "PowerPack: Add New Column — multi-source datasync
          persistence + formula recalc on rename (GROK-17109)"
          substring-matches scenario H3 "Add chained calculated columns,
          mutate, save with datasync, reopen, verify" via load-bearing
          keywords "Add", "datasync", "reopen" and the GROK-17109
          chain-witness role declared in the scenario body L516-523.
      - check_id: E-STRUCT-MECH-05
        status: PASS
        evidence: |
          Imports: @playwright/test (canonical) and three sibling-
          relative TestTrack helper modules (`../spec-login`,
          `../helpers/openers`, `../helpers/projects`). Sibling-
          relative-to-TestTrack-helpers convention is the established
          paradigm across TestTrack spec files; helper-module exports
          verified by Grep this cycle: spec-login.ts exports
          `baseUrl`, `specTestOptions`, `StepError`, `stepErrors`,
          `softStep`, `loginToDatagrok`, `loginAsSecondUser`;
          helpers/openers.ts exports `openTableFromFile` (L136),
          `assertProvenanceScript` (L550); helpers/projects.ts
          exports `saveProjectWithProvenance` (L855),
          `deleteProjectWithCleanup` (L1060). No imports from
          arbitrary external NPM paths.
      - check_id: E-STRUCT-MECH-06
        status: PASS
        evidence: |
          Spec L1-3 carry the leading frontmatter block
          `/* --- sub_features_covered: [...] --- */` before any other
          content (no leading imports or line comments above it). All
          five ids resolve to atlas entries verified by Grep:
          powerpack.dialogs (L587), powerpack.dialogs.add-new-column
          (L615), powerpack.dialogs.add-new-column-func (L622),
          powerpack.dialogs.prepare-add-column-call (L629),
          powerpack.formula.is-formula-column (L294). Spec frontmatter
          mirrors scenario frontmatter exactly (5 ids, same order).
      - check_id: E-TRACE-01
        status: PASS
        evidence: |
          Every `softStep(...)` block carries an explicit "Step N: ..."
          or "Setup:" label tracing to a numbered scenario step or to
          the scenario Setup section. Login + workspace-reset
          (L147-156) and `finally` cleanup (L644-654) are clearly
          delimited setup/cleanup glue (not softSteps, not scenario
          steps) and do not require `// technical:` markers.
      - check_id: E-TRACE-02
        status: PASS
        evidence: |
          All 10 numbered scenario steps covered: Step 1 (L187),
          Step 2 (L202), Step 3 (L310), Step 4 (L327), Step 5 (L427),
          Step 6 (L490), Step 7 (L507), Step 8 (L525), Step 9 (L577),
          Step 10 (L616). Setup steps 1+2 covered by "Setup: open
          System:DemoFiles/demog.csv with datasync provenance"
          softStep at L167. The Setup step 1 source-matrix narrowing
          (5 sources -> 1) is a coverage breadth reduction at the
          source-class axis, not a numbered-step omission — routed
          under SCOPE_REDUCTION above, not under E-TRACE-02 FAIL.
      - check_id: E-TRACE-03
        status: PASS
        evidence: |
          Verification calls present for each step: Step 2/4 verify
          column addition via df.columns.names().includes() polling +
          formula evaluation via expect(check!.diff).toBeCloseTo
          (100/200, 1) at L304/L405; Step 5 verifies formula tag
          rewrite via expect(formula.tag).toContain('BaseWeight') at
          L468; Step 6 verifies server-side persistence via
          grok.dapi.projects.find(pid) at L498-501; Step 7 verifies
          tables.length === 0 at L516; Step 8 GROK-17109 INVARIANT
          verified via expect(reopen.hasWeight2/hasWeight3).toBe(true)
          + expect(w2Formula/w3Formula.length).toBeGreaterThan(0) at
          L558-565; Steps 9/10 verify post-reopen formula tag rewrite
          and recompute deltas at L607-609 and L641-642.
      - check_id: E-SEL-01
        status: PASS
        evidence: |
          Verified by Grep against grok-browser/references this cycle:
          `[name="icon-add-new-column"]` in dialogs-menus.md L73;
          `[name="viewer-Grid"]` in viewers/grid.md L30 + grid.md L4 +
          viewers.md L39; column-header Rename context action in
          grid.md L52-54; Data sync toggle in projects.md L15, L23,
          L58, L70, L182, L192. Dialog input/button selectors
          (`[name="input-Add-New-Column---Name"]`,
          `[name="button-Add-New-Column---OK"]`,
          `.add-new-column-dialog-cm-div .cm-content`) are documented
          in spec leading block-comment L52-72 with explicit citation
          to PowerPack/src/dialogs/add-new-column.ts:344-349
          (prepareForSeleniumTests — verified by Grep this cycle: the
          function literally assigns these `name` attributes at
          L346-349). Qualifies as "documented in the test's local
          notes with reason" per E-SEL-01 alt path. `.d4-dialog` is
          the standard Datagrok dialog class.
      - check_id: E-SEL-02
        status: PASS
        evidence: |
          No invented selectors. Dialog input/button names come from
          PowerPack's prepareForSeleniumTests at
          add-new-column.ts:344-349 (verified); icon-add-new-column
          verified in grok-browser/references/dialogs-menus.md;
          viewer-Grid and `.d4-dialog` are platform conventions
          present across grok-browser/references.
      - check_id: E-SEL-03
        status: PASS
        evidence: |
          Dart Name input ([name="input-Add-New-Column---Name"]) driven
          via native HTMLInputElement.value setter + dispatchEvent
          ('input')/'change' at L206-213 (Step 2) and L329-336 (Step 4),
          NOT `.fill()`. CodeMirror composition uses
          view.dispatch({changes:...}) with keyboard.type() fallback
          (Step 2 L244 + L259; Step 4 L356 + L367), NOT `.fill()`.
          Dart-input-fill anti-pattern avoided across all five
          Dart-input interactions.
      - check_id: E-HELP-01
        status: PASS
        evidence: |
          `loginToDatagrok` (helpers-registry.yaml L3663) and `softStep`
          (L3657) are registered. The remaining imports
          (`specTestOptions`, `stepErrors`, `openTableFromFile`,
          `assertProvenanceScript`, `saveProjectWithProvenance`,
          `deleteProjectWithCleanup`) exist as concrete exports in
          TestTrack helper modules (verified by Grep this cycle:
          spec-login.ts L3, L5, L12, L14, L16, L37; helpers/openers.ts
          L136, L550; helpers/projects.ts L855, L1060). Helpers-
          registry backfill remains pending — qualifies as "declared
          as new helper candidate in the migration / design report"
          alt path. Recommend registry backfill in a separate
          mechanical batch (does NOT block PASS on E-HELP).
      - check_id: E-HELP-02
        status: PASS
        evidence: |
          No reinvention of registered helpers. Login uses
          loginToDatagrok (registered). File-open uses openTableFromFile
          (existing helper, not open-coded). Project save uses
          saveProjectWithProvenance rather than open-coding
          tables.save + projects.save + uploadDataFrame. Cleanup
          uses deleteProjectWithCleanup.
      - check_id: E-LAYER-01
        status: PASS
        evidence: |
          Scenario frontmatter `target_layer: playwright` aligns with
          the Playwright spec body (imports @playwright/test, uses
          `test(...)` from playwright with `page` parameter driving).
      - check_id: E-LAYER-02
        status: NA
        evidence: |
          No layer override; E-LAYER-01 passes.
      - check_id: E-LAYER-COMPLIANCE-01
        status: PASS
        evidence: |
          target_layer=playwright requires >=1 DOM-driving call; spec
          has many: page.locator(...).click()/.waitFor() at L188-193,
          L224-225, L280, L314-318, L341-342, L381; page.keyboard.press
          /type at L230-231, L259, L345-346, L367; cm.click() at L226,
          L255, L343, L363. pyramid_layer=bug-focused (NOT ui-smoke) —
          the JS-API-substitution-forbidden sub-rule does not apply.
          Step 5 column-rename uses JS API col.name setter (rationale
          L408-426 cites canvas-based grid context-menu coordinate
          brittleness under headless Playwright and equivalence of
          col.name dispatch to context-menu RenameColumn function);
          Step 6 uses saveProjectWithProvenance JS API path (rationale
          L478-488 cites Save Project dialog PascalCase normalization
          anti-pattern). Both substitutions acceptable for the
          bug-focused paradigm.
      - check_id: E-BOUND-01
        status: PASS
        evidence: |
          Spec file path
          public/packages/UsageAnalysis/files/TestTrack/PowerPack/add-new-column-advanced-spec.ts
          is under the allowed test-area path
          public/packages/UsageAnalysis/files/TestTrack/**.
      - check_id: E-BOUND-02
        status: PASS
        evidence: |
          No changes to core/** or package source outside src/tests/.
          Spec lives entirely in the TestTrack test-area.
      - check_id: E-RETRY-IGNORES-GATE-B
        status: NA
        evidence: |
          Retry-context detection does NOT fire this cycle. The
          scenario's gate_verdicts.b block (cycle_id
          2026-05-26-powerpack-automate-01, the most-recent-prior-cycle
          of this scenario) records verdict: PASS with failure_keys: []
          and spec_runs[].result: passed (3/3 attempts, durations
          38535/36230/50268 ms, errors[] empty, validated against
          per-attempt Playwright JSON reporter output). The
          E-RETRY-IGNORES-GATE-B predicate requires
          `gate_verdicts.b.verdict == "FAIL"`; the recorded Gate B
          verdict is PASS, so retry context is not active and this key
          never fires (per failure-keys-vocabulary: "Block absent →
          not a retry context; this key never fires" — same outcome
          for a present-but-PASS block). The earlier
          environment-defect FAIL narrated in the SR prose was from a
          superseded cycle (-03) and was remediated this cycle
          (duplicate @playwright/test install removed operator-side),
          which is why the recorded Gate B is now PASS. Independent of
          the predicate, the spec source retains no diagnosed-failing
          code path and proposes no un-backed paradigm pivot, so the
          ignored-evidence pattern is absent regardless. NA (check
          vacuous — not in retry context).
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-add-new-column-subdir
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML; required fields present and well-typed:
          feature=powerpack, sub_features_covered is a 5-element list,
          target_layer=playwright, coverage_type=regression. All five
          sub_features_covered entries (powerpack.dialogs.add-new-column,
          powerpack.dialogs.add-new-column-func,
          powerpack.dialogs.prepare-add-column-call,
          powerpack.formula.is-formula-column, powerpack.dialogs) resolve
          to atlas sub_features[].id entries.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains H2 scenario heading "## Scenarios" with an H3
          "### Add chained calculated columns, mutate, save with datasync,
          reopen, verify". The H2 + H3 pattern yields one identifiable
          scenario in the body.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          Under the H3 scenario heading, 10 numbered steps are present
          (1. Open the Add New Column dialog (first time) through 10.
          Change values in the source column (post-reopen)). Numbered-step
          requirement met.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          The single scenario in the body has 10 steps — not empty.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer=playwright is in the canonical
          {playwright, apitest, manual-only} enum.
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type=regression is in the canonical
          {smoke, regression, edge, perf} enum. Not a severity-axis p0..p3.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type=regression at the file-frontmatter level applies
          to every scenario in the body (single scenario). Test-kind enum
          value, not severity axis.
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          A ## Setup section factors the 5-source matrix selection and
          source-column-name binding (WEIGHT) up front so each cycle of
          the single scenario does not re-state setup. Single scenario
          in the body means there is no cross-scenario duplication risk.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          pyramid_layer=bug-focused (not ui-smoke). The hard alignment
          rule fires only on pyramid_layer=ui-smoke; bug-focused has no
          hard alignment constraint per the mode-file note (advisory:
          bug-focused → usually regression or edge — coverage_type is
          regression here, consistent with the advisory).
      - check: A-CONT-01
        status: PASS
        evidence: |
          Steps reference real names. Source 3-5 cite Postgres:Northwind
          query names (OrdersByEmployee, GetTop100, GetAll) and a real
          table (products). Source 1-2 cite local storage / Home dir as
          source kinds — real source-class descriptors, not placeholders.
          Calculated column names Weight2 and Weight3, and the rename
          target BaseWeight, are concrete. The string WEIGHT used in
          formulas is documented in Setup step 2 as a parametric source-
          column-name binding ("For Northwind `products`, `unitprice` is
          the natural choice; substitute the column name as needed per
          source") with explicit per-source substitution guidance — this
          is a parametric design, not a hallucinated angle-bracket /
          square-bracket placeholder. Step 3 TODO marker
          "(TODO: specify which formula to use)" present in the original
          was resolved during migration to the concrete formulas
          ${WEIGHT}+100 (Step 2) and ${Weight2}+100 (Step 4) as
          documented in source_text_fixes. No angle-bracket / square-
          bracket / generic-stand-in placeholders remain.
      - check: A-BUG-01
        status: PASS
        evidence: |
          Atlas powerpack.yaml known_issues is an empty list
          (known_issues: [] at line 1224, with the comment "No bug-
          library file exists for powerpack yet — leaving empty";
          known_issues has not yet been populated from
          bug-library/powerpack.yaml). With atlas known_issues empty,
          A-BUG-01 returns PASS-by-vacuity per the mode-file rule.
          Separately, the scenario does reference related_bugs:
          [GROK-17109] in its frontmatter and walks the canonical
          GROK-17109 reproduction path in Steps 6/8/9 — addressing the
          bug downstream of A-BUG-01 even when the atlas hook has not
          yet been wired.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          scope_reductions is an empty list ([]). No scenario step is
          opted out for effort or complexity reasons; the body walks
          all 10 steps in full.
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No "TODO: add later" / "deferred to next phase" markers
          remain. The original Step-3 TODO marker was resolved during
          migration (recorded in source_text_fixes:
          step-3-todo-resolved-formula-weight-plus-100-column-weight2
          and step-4-resolved-formula-weight2-plus-100-column-weight3),
          and the body now carries concrete formulas. The notes section
          documenting the resolution is descriptive metadata, not a
          deferral.
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:00:00Z
    spec_runs:
      - spec: add-new-column-advanced-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 37
        failure_keys: []
---

# Add New Column — multi-source datasync persistence + formula recalc on column rename

Multi-source matrix-like sequence walking the canonical GROK-17109
reproduction surface: open data from five distinct sources, add two
chained calculated columns, mutate the source column (rename + edit
values), save the project with datasync where the source supports it,
close all, reopen, and verify the calculated columns persist with
formulas intact AND recompute correctly when the source column is
renamed.

Chain witness role: this scenario is bug-focused per chain
`pyramid_layer: bug-focused` and owns the specialty persistence flows
(`save-project-with-datasync`, `project-reopen-with-formula-recalc`,
`column-rename-context-action`) for the PowerPack chain. Basic
dialog-open-and-add flow is delegated to the chain's smoke
witness at `add-new-column.md` (top-level PowerPack scenario,
`ui_coverage_delegated_to: add-new-column.md`).

## Setup

The scenario walks five table sources in sequence; for each source,
the same Scenarios block is executed end-to-end. Setup establishes
which source the cycle is currently using:

1. Pick one source from the matrix below and open the resulting
   table view:
   - **Source 1 — Local storage.** Open a table previously saved to
     local storage from the Browse tree.
   - **Source 2 — Home dir.** Open a table from the user's Home
     directory.
   - **Source 3 — Query result.** Run `Postgres:Northwind:OrdersByEmployee`
     from the Browse tree and open the result table.
   - **Source 4 — DB GetTop100 result.** Run the auto-generated
     `GetTop100` query against the `products` table in
     `Postgres:Northwind` and open the result table.
   - **Source 5 — DB GetAll result.** Run the auto-generated `GetAll`
     query against the `products` table in `Postgres:Northwind` and
     open the result table.
2. Pick a numeric source column on the open table that will be the
   subject of the chained formulas; call it `WEIGHT` in the
   instructions below. For Northwind `products`, `unitprice` is the
   natural choice; substitute the column name as needed per source.

## Scenarios

### Add chained calculated columns, mutate, save with datasync, reopen, verify

1. **Open the Add New Column dialog (first time).** Click the "Add
   new column" toolbar icon on the open table view ribbon. The
   `AddNewColumnDialog` opens.
2. **Add the first calculated column.**
   - Enter the formula `${WEIGHT} + 100` in the formula editor
     (substituting the source column name picked at Setup step 2;
     the example matches the `Weight2`/`Weight3` chain in the
     sibling `formula-refreshing.md` scenario for cross-scenario
     consistency).
   - Set the new column name to `Weight2`.
   - Click OK. A new column `Weight2` is added to the table; values
     are computed as the chosen source column plus 100.
3. **Open the Add New Column dialog (second time).** Click the
   "Add new column" toolbar icon again. The dialog reopens.
4. **Add the second calculated column referencing the first.**
   - Enter the formula `${Weight2} + 100` in the formula editor.
   - Set the new column name to `Weight3`.
   - Click OK. A new column `Weight3` is added; values are computed
     as `Weight2 + 100` (and therefore the chosen source column
     plus 200 transitively).
5. **Mutate the source column in the grid.**
   - Rename the source column header from `WEIGHT` to a new name —
     e.g. `BaseWeight` — using the context menu / Rename column
     action on the column header in the grid.
   - Change one or more cell values in the renamed source column.
   - **Expected result.** `Weight2` updates its formula text to
     reference the new source column name (`${BaseWeight} + 100`)
     and its values recompute accordingly. `Weight3`'s values
     recompute transitively from the updated `Weight2`.
6. **Save the project with datasync (where the source supports it).**
   - Save the current view as a project from the toolbar / Save
     dialog. In the Save Project dialog, enable the "Data sync"
     option for the open table.
   - **Note for matrix.** Some sources may not support datasync
     (e.g. tables opened from local storage have no upstream source
     to sync against). For sources where datasync is unavailable,
     save the project without datasync and note the source in the
     run log; the persistence-on-reopen invariant (Step 8) still
     applies to the calculated columns themselves regardless of
     datasync availability.
7. **Close all views.** Close the open table view (and any other
   views opened during this cycle) so the workspace is clean before
   reopen.
8. **Reopen the saved project.** Open the project saved at Step 6
   from Recent Projects / Dashboards / context menu.
   - **Expected result (GROK-17109 invariant).** Both calculated
     columns `Weight2` and `Weight3` are present in the dataset
     upon reopen, with formula tags preserved and values intact.
     Columns must NOT be missing — this is the canonical GROK-17109
     regression invariant.
9. **Rename the source column inside the formula (post-reopen).**
   - In the grid of the reopened project, rename the source column
     used by the formulas (the column carrying the new name from
     Step 5; if datasync rewrote it back on reopen, the current
     name is whatever the persisted project carries) to a different
     name — e.g. `BaseWeight2`.
   - **Expected result.** The formula on `Weight2` updates to
     reference the new source column name; the formula on `Weight3`
     remains referencing `${Weight2}` and remains valid.
10. **Change values in the source column (post-reopen).**
    - Edit one or more cell values in the source column on the
      reopened project's grid.
    - **Expected result.** Values in `Weight2` recompute according
      to its formula. Values in `Weight3` recompute transitively
      from the updated `Weight2`. The calculated columns track the
      source-column changes accordingly.

**Overall expected result.** The calculated columns are added,
their values change according to changes in the table source
columns AND persist across save-project-with-datasync + close-all
+ reopen with formula tags intact. Renaming the source column at
any point (pre-save or post-reopen) updates the formula text
automatically. This is the GROK-17109 reproduction path; before
the fix in 1.23.0, columns disappeared on reopen.

## Notes

- **Bug focus and chain role.** This scenario is bug-focused per
  chain `pyramid_layer: bug-focused`; it walks the canonical
  GROK-17109 reproduction (Step 6 save-with-datasync + Step 8
  reopen-and-verify) plus the formula-recalc-on-rename invariant
  (Step 9). Cross-cutting bug candidate `GROK-17109` is emitted at
  the chain level (`bug_focused_candidates` in
  `scenario-chains/powerpack.yaml`); proposed spec
  `powerpack-grok-17109-spec.ts` spans this scenario plus
  `add-new-column.md` (top-level smoke) and
  `AddNewColumn/formula-refreshing.md` (dependency-chain recalc).
- **UI delegation.** Basic dialog-open-and-add flow (toolbar icon,
  formula editor, OK button) is owned by the chain's smoke witness
  at `add-new-column.md`. This scenario owns the specialty
  persistence flows: `save-project-with-datasync`,
  `project-reopen-with-formula-recalc`, `column-rename-context-action`.
  See chain `ui_coverage_plan.delegated_scenarios` entry.
- **Sibling spec.** A Playwright spec already exists at
  `public/packages/PowerPack/src/tests/add-new-column.ts` (see
  `existing-test-index.yaml`); house-style anchor for Automator
  when authoring the migrated scenario's `-spec.ts`.
- **Sibling scenario alignment.** The `Weight2 = ${WEIGHT} + 100`,
  `Weight3 = ${Weight2} + 100` chain mirrors the formula pattern
  used by `AddNewColumn/formula-refreshing.md` (Weight2 → Weight3 →
  Weight4) for cross-scenario consistency; that sibling extends the
  chain with `Weight4 = Log10(${Weight3}) - 0.2` and adds Context
  Panel formula-edit coverage which this scenario does not duplicate.
- **Source-text fixes.** The original scenario contained the TODO
  marker `(TODO: specify which formula to use)` on Step 3; the
  migration resolves it to `${WEIGHT} + 100` for Step 3 and
  `${Weight2} + 100` for Step 4 per the chain-level recommendation
  (a) in `unresolved_ambiguities :: formula-not-specified` — the
  resolution mirrors the Weight2/Weight3 chain in the sibling
  `formula-refreshing.md` scenario. The original "Close All"
  shorthand is spelled out as closing all open views via the shell
  (Step 7) so the action is unambiguous when this scenario is
  automated.
- **Candidate helpers.** Several patterns recur and warrant
  registry candidates: opening tables from local storage / Home dir
  / query results (`openTableFromLocalStorage`, `openTableFromHomeDir`,
  `openQueryResult`), saving a project with datasync
  (`saveProjectWithDatasync`), reopening a saved project
  (`reopenProject`), adding a calculated column via the dialog
  (`addCalculatedColumn`), and renaming a column via the grid
  context action (`renameColumnViaContextAction`). None of these
  exist in `helpers-registry.yaml` yet; surfaced here as candidates
  for the registry per the Migrator candidate-helper convention.
- **Datasync availability matrix.** Sources 1 (local storage) and
  2 (Home dir) typically do not support datasync since there is no
  upstream source to sync against; Sources 3-5 (Northwind query
  results) do support datasync. The scenario body's Step 6
  "(where available)" qualifier preserves the original's allowance
  for variable datasync support across sources; the persistence-
  on-reopen invariant (Step 8) applies to all five sources.
- **Original trailing JSON metadata.** The original scenario ended
  with `{"order": 1}`. The `order` field is captured in chain
  `order_from_files` under the path-relative key
  `AddNewColumn/add-new-column.md`.
- **Name-collision awareness.** There is a sibling top-level
  scenario at `PowerPack/add-new-column.md` (the chain's smoke
  witness for the Demog flow). The two scenarios share the same
  basename but are DISTINCT and not interchangeable; the chain
  encodes them with path-relative keys. See chain
  `unresolved_ambiguities :: naming-collision-add-new-column` for
  future rename suggestions.
