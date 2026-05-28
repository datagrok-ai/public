---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-functions-panel
  - add-new-column-functions-sort-by-name
ui_coverage_delegated_to: add-new-column.md
ui_companion: functions-sorting-ui.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/functions-sorting.md
migration_date: 2026-05-20
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
realized_as:
  - functions-sorting-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-functions-sorting
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:05:00Z
    failure_keys: []
    scope_reduction_proposal: |
      Independent fresh-context Critic E re-evaluation under
      inputs.cycle_id = 2026-05-28-powerpack-automate-02. Same-cycle
      retry context per gate_verdicts.b.cycle_id == inputs.cycle_id
      (Validator FAIL at 2026-05-28T10:10:00Z, failure_keys
      [B-RUN-PASS, B-STAB-01], 3/3 attempts, 26s). First OR branch
      of mode-file §"Gate B awareness" detection rule satisfied.
      Full spec body and full scenario re-read; supporting context
      (helpers-registry, atlas slice, grok-browser references,
      decision-log failed_attempts) re-verified independently. Prior
      gate_verdicts.e (cycle-01) NOT consulted for this verdict;
      reached the same SCOPE_REDUCTION place by independent reasoning.

      Mechanical fail-fast (E-STRUCT-MECH-01..06) — all PASS:
        01: spec exists at functions-sorting-spec.ts (same dir).
        02: TS parses; braces balanced; single test(...) opens at
          L622, closes at L1413; imports well-formed.
        03: exactly one test(...) block at L622.
        04: test name "PowerPack: Add new column functions-panel
          sorting (SPGI — by type, by name, sticky)" contains
          substring anchors "sort"/"by name"/"sticky" from the
          scenario heading. Substring traceability satisfied.
        05: imports `@playwright/test` (canonical) + relative
          `../spec-login`. softStep (helpers-registry L3657) and
          loginToDatagrok (L3663) registered; specTestOptions
          (spec-login.ts L5) and stepErrors (L14) are sibling-
          module locals — house style. No non-canonical import.
        06: leading `/* --- sub_features_covered: [...] --- */`
          block at L1-3 (no leading imports / comments above).
          Both `powerpack.dialogs.add-new-column` (atlas L608) and
          `powerpack.dialogs.add-new-column-func` (atlas L615)
          resolve to feature-atlas/powerpack.yaml ids (Grep-verified).

      Retry-context predicate (mode-file §"Gate B awareness"):
        gate_verdicts.b.verdict == "FAIL" AND gate_verdicts.b.cycle_id
        == "2026-05-28-powerpack-automate-02" == current inputs.cycle_id.
        Predicate satisfied — same-cycle retry, first OR branch.

      E-RETRY-IGNORES-GATE-B — NOT fired:
        Validator's diagnosis (gate_verdicts.b L411-420):
        [B-RUN-PASS, B-STAB-01], 3 attempts, 26s. The B block carries
        no spec_runs[].flake_evidence pinpointing the failing
        assertion, so the cross-check is against the failure-key
        semantics (B-RUN-PASS = ≥1 assertion executed and failed;
        26s rules out a B-STAB-04 timeout surface). Cross-check of
        the current spec body against this diagnosis:
          (i) The current cycle (2026-05-28) DID author a substantive,
              identifiable diff — NOT cosmetic. L1211-1226 records a
              fresh dev.datagrok.ai 2026-05-28 MCP recon re-confirming
              the canvas-click trigger produces FOUR distinct
              function-list orderings (numeric / default / boolean /
              Molecule families) and re-establishing the row→family
              map (rows 0/12 numeric, 1-6/13 default, 7 boolean, 8-11
              Molecule). The earlier "canvas untestable" claim is
              explicitly retracted as stale (L1224-1226, L914-919).
          (ii) Step 5 settle logic rewritten (L1316-1328): the
               change-detector keyed on top-2 entries is replaced by
               a poll-for-two-consecutive-identical-reads settle. The
               inline rationale (L1310-1317) cites the exact B-RUN-PASS
               hazard: numeric-input top-2 (Abs, Acos) is IDENTICAL to
               alphabetical top-2 (Abs, Acos), so a top-2 change-
               detector returns a mid-render snapshot, and Step 6 then
               diffs Step 6's byte-for-byte sticky baseline against a
               non-settled list — a deterministic assertion failure.
               Step 6 (L1362-1386) now diffs against the settled Step 5
               baseline. This is a direct, MCP-backed fix to a
               plausible B-RUN-PASS surface for [B-RUN-PASS, B-STAB-01].
          (iii) The Steps 3/4 probe approach (findRowProducingDistinctOrder
                + clickColumnRowByIdx) and the retry-1 setup-phase
                semType-race poll are preserved — those addressed prior
                cycles' distinct failure surfaces and are not regressed.
        The change is an identifiable diff against the prior failing
        path, MCP-empirically backed, targeting a plausible B-RUN-PASS
        surface; it is NOT cosmetic and NOT a paradigm pivot lacking
        backing (canvas-click trigger retained; the JS-API selectColumn
        path is explicitly dead/unused at L905-912, L988-993). Per
        mode-file §"Gate B awareness" final paragraph: partial /
        good-faith retries route via PASS / SR / EG; Critic E does NOT
        enforce convergence (that is Validator's job at Gate B on the
        next run). No do_not_retry directive in decision-log
        failed_attempts forbids the canvas-click approach for
        functions-sorting (the do_not_retry entries are all `legend`).

      Content checks:
        - E-TRACE-01: PASS. softStep bodies for Steps 2/3/4/5/6 trace
          1:1 to scenario numbered steps; the Step 1 sanity block
          (L715-746) is the scenario Step 1 Verify (semType / column-
          name reads). The pre-Step-2 onDialogShown subscription
          (L754-820) is scaffolding for the now-dead selectColumn
          helper — documented inline as retained-but-unused (L905-912),
          not orphaned untracked code.
        - E-TRACE-02: PASS at the abstract-flow level (SR-02 loosening
          on Steps 3/4 example-specificity). All 6 numbered scenario
          steps have ≥1 corresponding test step.
        - E-TRACE-03: PASS at the abstract-assertion level (SR-02
          loosening). Step 3/4 verifies map to order-changed assertions
          (L1249, L1265-1266); specific-family-on-top demoted to
          log-only audit (L1252-1253, L1267-1268). Steps 1/2/5/6
          assertions at full scope.
        - E-SEL-01/02: PASS. `[name="icon-add-new-column"]` at
          grok-browser/references/dialogs-menus.md:73 (Grep-verified).
          `.d4-dialog`, `.d4-menu-popup`, `.d4-menu-item` in
          grok-browser widgets/dialog.md + widgets/menu.md (Grep-
          verified). `[name="icon-sort-alt"]` is NOT in grok-browser
          references but is cited verbatim by the scenario body Step 5
          (scenario L491) AND documented in the spec header (L51-53)
          with a Dart source citation (functions_view.dart:351,
          `Icons.faSolid('sort-alt', ...)`) — satisfies the
          local-notes branch of E-SEL-01; not an invented selector.
          Dialog-internal selectors are documented with source-file
          line citations in spec header L34-65.
        - E-SEL-03: PASS. No Dart-input `.fill()`. Triggers are
          Playwright `.click()` on real DOM locators or synthetic
          canvas.dispatchEvent(MouseEvent) inside page.evaluate —
          canvas-DOM driving, not `.fill()`-style Dart-input typing.
        - E-HELP-01/02: PASS. softStep, loginToDatagrok registered;
          specTestOptions, stepErrors sibling-module locals. Inline
          helpers (readFunctionOrder, clickColumnRowByIdx,
          findRowProducingDistinctOrder, waitForOrderChange,
          pickColumnsBySemType, dead selectColumn) are spec-local; no
          reinvention of registered helpers.
        - E-LAYER-01/02: PASS. target_layer=playwright matches the
          Playwright body shape; no override declared; decision-log
          layer_decisions: [] (vacuous).
        - E-BOUND-01/02: PASS. Spec under TestTrack/PowerPack/ — an
          allowed test-area path. No core/** or non-test package
          source edits.

      E-LAYER-COMPLIANCE-01:
        - target_layer=playwright requires ≥1 DOM-driving call:
          abundantly satisfied (page.locator/.click/.waitFor,
          dlg.locator across Steps 2/5/cleanup; canvas dispatchEvent
          inside page.evaluate against .add-new-column-columns-grid).
        - pyramid_layer=ui-smoke sub-rule: JS-API substitution
          forbidden for owned UI flows. The sort-by-type trigger is
          canvas synthetic-MouseEvent dispatch on a real DOM canvas
          element — DOM-event-chain driving, NOT the JS-API
          substitution (`columnsDf.currentRowIdx = N`) the rule
          forbids; that JS-API path is the explicitly-dead selectColumn
          helper (L988-993). sort-by-name (Steps 5/6) is pure Playwright
          Locator clicks. Strict regex-based mechanical reading:
          satisfied. PASS.

      SCOPE_REDUCTION basis (load-bearing): SR-02 (assertion-family-
      specificity loosening) is REFLECTED in the spec body (Steps 3/4
      assert order-changed-only with specific-family-on-top demoted to
      log-only) but is NOT PERSISTED into the scenario frontmatter's
      scope_reductions[] list (scenario L20: `scope_reductions: []`).
      Mode-file §"Verdict guidance" for SCOPE_REDUCTION: "if the spec
      covers fewer assertions than the scenario explicitly because of a
      documented technical limitation" — this is exactly that case
      (platform-catalogue drift between scenario-body example function
      families and the platform's current top-of-list output). The
      reduction is fully reasoned in the spec header and MCP-backed,
      but the scenario-frontmatter declaration is the durable record
      consumed downstream.

      Proposed reductions (to be persisted by the Migrator / operator
      into scenario frontmatter scope_reductions[]):

        scope_reductions:
          - id: SR-01
            check: ui-smoke-trigger-canvas-affordance
            rationale: |
              The owned UI flow `add-new-column-functions-sort-by-type`
              is triggered by canvas synthetic-MouseEvent dispatch on
              the canvas inside .add-new-column-columns-grid rather than
              a Playwright Locator `.click()`. The columns widget is a
              canvas-based DG.ColumnGrid.popup (widgetMode: true;
              PowerPack/src/dialogs/add-new-column.ts:1079) — there is
              no per-row DOM element a Locator can target. Live MCP
              recon (2026-05-26 / 2026-05-27 / 2026-05-28) confirmed
              (a) JS-API columnsDf cannot resolve for the popup-mode
              ColumnGrid from a Playwright spec, and (b) canvas
              synthetic-MouseEvent dispatch DOES trigger the
              functions-list re-sort. sort-by-name (Step 5) and the
              sticky-sort popup mechanic (Step 6) remain end-to-end
              DOM-driven via Playwright Locator. SR-01 is informational
              under the strict reading of E-LAYER-COMPLIANCE-01 (canvas
              dispatchEvent is DOM-driving, not JS-API substitution);
              persistence makes the canvas-affordance technical
              constraint discoverable for future refresh / atlas
              curation.
            verdict_status: SCOPE_REDUCTION
          - id: SR-02
            check: ui-smoke-assertion-function-family-specificity
            rationale: |
              Scenario Steps 3 and 4 cite specific Molecule-input and
              numeric-input function names as the expected top-of-list
              after the column-selection trigger (canonicalize,
              convertMolNotation, getCLogP, getDescriptors for Structure;
              Abs, Acos, Asin, Atan, Atan2 for numeric). Live MCP recon
              (dev.datagrok.ai 2026-05-26 / 2026-05-28) shows the
              platform's current function catalogue does NOT place these
              examples in the top-10 after the corresponding column
              click (after Structure click, top-10 is [BDE_prediction,
              Column Exists, Contains, createProperty, createTemplate,
              ...]; canonicalize is at position 347 of 455). The
              canvas-row → source-df-column mapping is also non-linear
              (rows 8-11 produce the Molecule family, not row 1).
              Steps 3/4 assertions are loosened to order-changed-only
              via findRowProducingDistinctOrder; specific-function-
              family-on-top is a log-only audit signal. Sort-by-name
              (Step 5) and sticky-sort (Step 6) assertions remain at
              full scope — those checks are universal and pass against
              the platform's actual behavior. The scenario-cited
              function-family examples remain useful as documentation /
              atlas-curation candidates for future scenario refresh.
            verdict_status: SCOPE_REDUCTION

      Note: SR-02 is the load-bearing reduction (assertion loosening
      from specific-family-on-top to order-changed-only). SR-01 is
      informational under the strict E-LAYER-COMPLIANCE-01 reading
      (canvas dispatchEvent is DOM-driving) and is not a hard
      prerequisite for that check to PASS.
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-functions-sorting
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
    review_round: 1
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
        evidence: |
          Frontmatter parses as YAML and contains all four required fields:
          feature=powerpack, sub_features_covered=[powerpack.dialogs.add-new-column,
          powerpack.dialogs.add-new-column-func], target_layer=playwright,
          coverage_type=smoke.
      - check: A-STRUCT-MECH-02
        status: PASS
        evidence: |
          Body contains a scenario heading at level 3 ("### Functions panel
          re-sorts by column-input type, then by name, then sticks") under
          the "## Scenarios" section.
      - check: A-STRUCT-MECH-03
        status: PASS
        evidence: |
          The single scenario contains six numbered steps (1. through 6.)
          covering dataset open, dialog open, type-sort by chem column,
          type-sort by numeric column, switch to By-name sort, and the
          sticky-sort contract.
      - check: A-STRUCT-MECH-04
        status: PASS
        evidence: |
          The scenario is non-empty — six numbered steps present with
          Verify sub-clauses for each.
      - check: A-STRUCT-MECH-05
        status: PASS
        evidence: |
          target_layer=playwright is in the allowed enum {playwright,
          apitest, manual-only}.
      - check: A-STRUCT-MECH-06
        status: PASS
        evidence: |
          coverage_type=smoke is in the allowed enum {smoke, regression,
          edge, perf}.
      - check: A-STRUCT-03
        status: PASS
        evidence: |
          coverage_type=smoke is set at file-frontmatter level and applies
          to the single scenario in the body. No severity-axis (p0..p3)
          value used.
      - check: A-STRUCT-04
        status: PASS
        evidence: |
          Setup (clean Datagrok session, depends_on=[] per chain) is
          factored into the "## Setup" section. With only one scenario in
          the file there is no cross-scenario duplication to flag.
      - check: A-LAYER-ALIGN-01
        status: PASS
        evidence: |
          Frontmatter pyramid_layer=ui-smoke is paired with coverage_type=smoke
          per the alignment rule. Chain YAML powerpack.yaml line 322 confirms
          pyramid_layer=ui-smoke for AddNewColumn/functions-sorting.md.
      - check: A-CONT-01
        status: PASS
        evidence: |
          Scenario references real column names (Structure, Chemical Space X,
          Chemist), real file path (System:DemoFiles/chem/SPGI.csv), real
          selector ([name="icon-sort-alt"]), and real function names
          (canonicalize, convertMolNotation, getCLogP, Abs, Acos, Atan2).
          No angle-bracket / square-bracket placeholders or generic
          stand-ins detected.
      - check: A-BUG-01
        status: PASS
        evidence: |
          PASS-by-vacuity. Atlas powerpack.yaml known_issues field is
          empty (line 1224: known_issues: []), so no needed-coverage bugs
          exist for this gate to check. related_bugs=[] in the scenario
          frontmatter is consistent with this.
      - check: A-MERIT-01
        status: PASS
        evidence: |
          No scenario opted out for effort or complexity. The single
          scenario is fully present with 6 numbered steps.
      - check: A-MERIT-02
        status: PASS
        evidence: |
          No unprompted deferrals. The Notes section references a
          curator candidate for next atlas regen (functions-panel sort-icon
          menu not currently captured as atlas interactions[] entries) but
          this is informational about atlas curation, not a deferral of
          scenario work.
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T10:48:30Z
    spec_runs:
      - spec: functions-sorting-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 30
        failure_keys: []
---

# Add New Column dialog — functions-panel sorting

Focused UI scenario verifying the functions-panel sorting behavior in
the Add New Column dialog: (a) selecting a column re-sorts the functions
list so functions whose input parameter type matches the selected
column appear at the top (Molecule-input functions for chem columns,
numeric-input functions for numeric columns); (b) toggling the sort
icon to "By name" re-sorts the list alphabetically; (c) once "By name"
is active, clicking different columns does NOT change the function order.

`pyramid_layer: ui-smoke` per `scenario-chains/powerpack.yaml` rev 1
(Rule 1 fallback — single short UI flow over one source `SPGI`,
`simple` classification, no related bugs). NOT the chain's elected
smoke witness; `coverage_type: regression` is used here (the chain
reserves `coverage_type: smoke` for the elected smoke scenario at
top-level `add-new-column.md`, which exercises the dialog UI
holistically).

`ui_coverage_delegated_to: add-new-column.md` — this scenario owns
the functions-panel sort-by-type / sort-by-name specialty flows on
its responsibility list and delegates the basic dialog-open-and-add
flow to the smoke witness.

`related_bugs: []` — no curated bug in
`bug-library/powerpack.yaml` intersects this scenario's
`sub_features_covered` for the functions-panel-sorting surface.
(GROK-17109 and GROK-17004 affect the same `powerpack.dialogs.add-new-column`
sub-feature but their reproduction surfaces — formula-recalc-on-rename
across save-with-datasync, and complex-paste-handler crash — are not
exercised by the functions-list sort behaviour. See chain rev 1
`bug_focused_candidates[]` for the canonical cross-cutting spec
proposals.)

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining (`depends_on: []` per chain rev 1).

## Scenarios

### Functions panel re-sorts by column-input type, then by name, then sticks

1. Open the SPGI dataset (`System:DemoFiles/chem/SPGI.csv`) — for
   example via **Browse** > **Files**, **File** > **Open**, or the
   equivalent JS-API loader. **Verify:** the SPGI grid renders with
   its chem `Structure` column (semType `Molecule`) plus numeric
   and string columns (e.g. `Chemical Space X`, `Chemist`).
2. Open the **Add New Column** dialog — click the `Add New Column`
   icon on the table view toolbar (or **Edit** > **Add New Column**
   from the top menu). **Verify:** the Add New Column dialog opens
   with its formula editor (CodeMirror), columns list, functions
   panel, and preview grid visible. The functions panel starts in
   its default sort mode ("By relevance").
3. Click the `Structure` column in the dialog's columns list.
   **Verify:** the functions list on the right re-sorts so that
   functions with `Molecule`-type input parameters appear at the
   top (for example, `canonicalize(molecule)`,
   `convertMolNotation(molecule, ...)`, `convertMoleculeNotation(molecule, ...)`,
   `getCLogP(smiles)`, `getDescriptors(molecules, ...)` precede
   functions whose first input is numeric or string).
4. Click a numeric column (e.g. `Chemical Space X`, type `double`)
   in the columns list. **Verify:** the functions list re-sorts so
   functions whose first input parameter is numeric appear at the
   top (for example, `Abs(x)`, `Acos(x)`, `Asin(x)`, `Atan(x)`,
   `Atan2(a, b)` precede chem-input or string-input families).
   Repeat for any other column type as available and **verify:**
   the matching-parameter family is on top in each case.
5. Click the sort icon on the top right corner of the functions
   panel (`[name="icon-sort-alt"]` — the two blue arrows). **Verify:**
   a popup menu appears with options "By name" and "By relevance".
   Select "By name". **Verify:** the functions list re-sorts
   alphabetically (for example, `Abs`, `Acos`, `Add`, `And`, `Asin`,
   `Atan`, `Atan2`, `Avg`, `BinByDateTime`, `BinBySpecificLimits`,
   `Boolean`, `Call`, ... appear in alphabetical order from the top).
6. With "By name" sort still active, click several different columns
   in the columns list in succession (for example `Structure`, then
   `Chemical Space X`, then `Chemist`). **Verify:** the alphabetical
   function order does NOT change — clicking columns no longer
   re-orders the functions list once "By name" sort is selected
   (sticky-sort contract).

## Notes

- **target_layer: playwright** — chosen because a sibling
  `functions-sorting-spec.ts` already exists at the playwright layer
  (per `existing-test-index.yaml` line 33452, layer `playwright`).
  All other Add-New-Column siblings (`add-new-column-spec.ts`,
  `autocomplete-spec.ts`, `hints-spec.ts`, `highlight-spec.ts`,
  `input-functions-spec.ts`, `formula-refreshing-spec.ts`) are also
  Playwright specs, so the house style is consistent.
- **coverage_type: regression** — per the chain's smoke-scenario
  election (`ui_coverage_plan.smoke_scenario: add-new-column.md`),
  the top-level `add-new-column.md` owns the `smoke` slot. This
  narrower per-flow scenario lands in `regression`. (Per the
  migration prompt's `coverage_type` heuristic: `smoke` for one
  fast happy-path test per feature; `regression` is the default
  for the rest.)
- **Step decomposition note.** The original scenario's Step 3 packs
  two sub-cases into one numbered step: (a) Structure column →
  Molecule-input functions on top; (b) change to a non-Structure
  column and verify the matching-input family is on top. The
  migrated body splits these into Step 3 (Structure / chem) and
  Step 4 (numeric column). Both sub-cases of the original Step 3
  are preserved verbatim — no silent drop. The chain's
  `ui_coverage_responsibility:` for this scenario already
  enumerates both `add-new-column-functions-sort-by-type` and
  `add-new-column-functions-sort-by-name`, so the split aligns
  with the chain's coverage intent.
- **Helpers (already in registry, available for downstream
  Automator):** `softStep`, `loginToDatagrok`, `specTestOptions`,
  `stepErrors` from
  `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts` —
  used by the existing `functions-sorting-spec.ts` per the sibling
  test index entry's `helpers_called: [spec-login]`.
- **Bug-library status:** consulted —
  `bug-library/powerpack.yaml` exists and was scanned; no curated
  bug intersects this scenario's `sub_features_covered` for the
  functions-panel sort surface. Cross-cutting candidates GROK-17109
  and GROK-17004 are emitted at the chain level
  (`bug_focused_candidates[]`) but their `spans` reference
  `AddNewColumn/add-new-column.md`, `AddNewColumn/formula-refreshing.md`,
  and `AddNewColumn/highlight.md` — NOT this scenario.
- **Decision log status:** queried — no
  `failed_attempts WHERE feature == powerpack` entries that touch
  the `AddNewColumn/functions-sorting.md` migration. No retry-skip
  approaches apply.
- **Atlas linkage (`derived_from:` provenance):**
  - `powerpack.dialogs.add-new-column-func` interactions "click
    Add New Column toolbar icon → opens AddNewColumnDialog" and
    "Edit | Add New Column top-menu → opens AddNewColumnDialog"
    are code-derived from
    `public/packages/PowerPack/src/package.ts#L405`.
  - `powerpack.dialogs.add-new-column` is the dialog class with
    the columns list + functions panel + preview grid; source
    anchor `public/packages/PowerPack/src/dialogs/add-new-column.ts#L98`.
    The functions-panel sort-by-type behaviour (matching first
    input parameter type to the selected column) and the sort-icon
    menu ("By name" / "By relevance") are surfaces of this class
    not currently captured as atlas `interactions[]` entries —
    surfaced as a curator candidate for the next atlas regen, but
    not blocking for this migration.
- **Edge cases.** The original scenario implies an edge case in
  Step 5 ("Try to click on columns — function order should not
  change"): the alphabetical-sort sticky-contract. Preserved as
  Step 6 of the migrated body. Step 3's "Change column to some
  other and check that function are sorted by matching parameters"
  generalises across column types; the migrated Step 4 picks a
  representative numeric column (`Chemical Space X`) and notes the
  pattern applies to any column type — matches the
  `functions-sorting-run.md` evidence, which records `string`
  (`Chemist`) as covered by the same assertion.

---
{
  "order": 5
}
