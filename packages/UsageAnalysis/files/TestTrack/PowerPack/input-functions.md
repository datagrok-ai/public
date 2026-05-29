---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-function-plus-icon
ui_coverage_delegated_to: add-new-column.md
ui_coverage_split_to:
  - input-functions-ui.md
realized_as:
  - input-functions-spec.ts
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/input-functions.md
migration_date: 2026-05-20
source_text_fixes:
  - add-new-column-dialog-title-case-normalization
  - text-filed-to-text-field-typo
  - fucntion-typo-corrected-to-function
  - drag-n-drop-phrase-normalization
candidate_helpers: []
unresolved_ambiguities:
  - sort-icon-glyph-described-not-named
  - id-column-type-classification
scope_reductions:
  - id: SR-01
    check: ui-smoke-drag-drop-affordance
    rationale: |
      Owned ui-smoke flow `add-new-column-function-drag-drop` exercised via
      the insertIntoCodeMirror END STATE
      (`document.execCommand('insertText', '<name(arg)>')` on the focused
      CM6 .cm-content), identical to what the platform DnD produces. MCP
      recon 2026-05-28 on dev.datagrok.ai confirmed the drag-drop UI leg is
      a genuine affordance gap: function rows are NOT HTML5-`draggable`
      (`d4-link-label` / `TR.d4-current-object`, draggable=false), Datagrok
      uses Dart pointer-event DnD via the `_dndContext` registry seeded
      source-side, so Playwright native dragTo (HTML5 drag) does not fire it
      and synthetic DragEvent corrupts the editor (inserts the raw
      text/plain name). The plus-icon `+` insertion leg IS genuine
      DOM-driving (trusted .click() on [name="icon-plus"] propagates).
      `usedFallback: true` surfaces the affordance gap via console.warn for
      reviewer.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: ui-smoke-column-grid-selection-column-name-relaxation
    rationale: |
      REVISED from the cycle 2026-05-26 round-2 SR-02
      (column-grid-selection-affordance). MCP recon 2026-05-28 on
      dev.datagrok.ai REFUTED the round-2 conclusion that column-grid
      selection is "genuinely NOT scriptable": a synthetic-MouseEvent triple
      (mousedown+mouseup+click) dispatched on the column-grid's LAST canvas
      (overlay) DOES select a column and DOES fire the dialog's
      `columnsDf.onCurrentRowChanged` (add-new-column.ts:1095), which sets
      `selectedColumn` (auto-bind) AND `sortByColType` (re-sort). The
      positive auto-bind contract is therefore RESTORED — the spec asserts
      `Chem:getCLogP(${<MoleculeCol>})` and `Abs(${<NumericCol>})` for real
      (no console.warn downgrade). The ONLY relaxation: the canvas row ->
      column mapping is non-linear (the popup groups columns by inferred
      input family) and the SPECIFIC `Structure` column (source idx 1) is
      NOT among the 14 visible canvas rows on SPGI (wheel-scroll did not
      reach it; "Search column" filter + row-0 click did not re-map the
      hit-test). So the auto-bind type-match invariant is asserted against
      whichever Molecule / numeric column the canvas PROBE selects
      (Core/R1/R3/R100 / Chemical Space X|Y on SPGI), not literally
      `${Structure}`. The scenario's INTENT (type-match auto-bind for a
      Molecule / numeric column) is fully verified; only the literal example
      column name relaxes. The Step 10 negative contract (no auto-bind on
      type mismatch -> `Abs(num)`) is asserted against a string-column probe
      (the SPGI `Id` column is string/semType=null), and holds.
    verdict_status: SCOPE_REDUCTION
related_bugs: []
gate_verdicts:
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T16:05:00Z
    failure_keys: []
    scope_reduction_proposal: |
      Fresh independent Gate E review of input-functions-spec.ts.

      E-STRUCT-MECH-01..06 PASS: spec present at
      PowerPack/input-functions-spec.ts; balanced braces, well-formed
      test() signature (line 463), imports for all used symbols; one
      test() block (line 463); test name "Add new column — function
      insertion (plus icon, drag-and-drop, auto-bound column parameter)"
      substring-overlaps the scenario heading "Insert functions via +
      icon and drag-and-drop; auto-bind columns by type match"; imports
      only @playwright/test + ../spec-login (registered sibling helper
      module); leading /* --- sub_features_covered: [...] --- */ block at
      lines 1-3 before any other content, non-empty, both ids resolve to
      atlas sub_features (powerpack.yaml lines 608 and 615).

      E-TRACE-01..03 PASS: all 10 numbered scenario steps map to softStep
      blocks (Steps 2-10 + 7a/7b + 10b); Step 1 is the SPGI-open
      precondition asserting cols.contains('Structure') + semType Molecule
      (lines 509-510); each Verify maps to an expect(...) JS-API assertion
      (col/semType checks, formula regexes, preview presence); probe/clear/
      read helpers are documented technical glue for the non-linear canvas
      row->column mapping, not orphaned code. E-SEL-01..03 PASS: no
      grok-browser/references/powerpack.md exists (confirmed via glob), so
      E-SEL-01's local-notes-with-reason clause applies — every selector
      is documented in the spec header (lines 54-73) with add-new-column.ts
      source citations; no invented selectors; no .fill() on Dart inputs.
      E-HELP-01/02 PASS: softStep + loginToDatagrok are registered
      (helpers-registry.yaml lines 3735, 3741); specTestOptions/stepErrors
      are config/error-array exports of the same spec-login.ts module
      (lines 5, 14), not reinvented helpers. E-BOUND-01/02 PASS (spec under
      TestTrack/PowerPack/, no core/** touched).

      E-LAYER-01/02 + E-LAYER-COMPLIANCE-01 PASS: target_layer playwright
      with many DOM-driving calls (page.locator/.click/.hover/.evaluate).
      Pyramid sub-rule (pyramid_layer ui-smoke): scenario
      ui_coverage_responsibility lists ONE owned flow
      (add-new-column-function-plus-icon), exercised via a genuine TRUSTED
      .click() on [name="icon-plus"] (clickPlusIcon, lines 289-314) — no
      JS-API substitution. Calls(>=1) >= flows(1). PASS.

      RETRY CONTEXT: gate_verdicts.b.verdict == FAIL, cycle_id
      2026-05-28-powerpack-automate-02 == current inputs.cycle_id,
      failure_keys [B-RUN-PASS, B-STAB-01] — same-cycle retry, so
      E-RETRY-IGNORES-GATE-B was evaluated. It does NOT fire. The new spec
      makes two real, identifiable, non-cosmetic diffs in the exact failing
      region: (a) adds VALIDATION_TYPES_MAPPING (lines 199-207) applied in
      dragFunctionOntoEditor's auto-bind position search (line 429,
      mapped.includes(args.sel.type)), so double/int/qnum columns bridge to
      Abs's `num` input — addressing the cited B-RUN-PASS assertion miss;
      (b) moves editor read/clear/insert from the desync-prone
      document.execCommand contenteditable path to the live CM6 EditorView
      dispatch (readEditorDoc line 228, clearEditor line 252,
      dragFunctionOntoEditor line 442), addressing the text-accumulation
      symptom. These are good-faith changes to the failing path (not a
      retained failing path, not cosmetic-only). Whether the fix holds at
      runtime is Validator's job at Gate B next run; Critic E does not
      enforce convergence. (Decision log records a prior
      E-RETRY-IGNORES-GATE-B FAIL in cycle 2026-05-26-powerpack-automate-02
      for this feature — that earlier round genuinely ignored the
      diagnosis; this round does not.)

      Two acceptable scope reductions are recorded on the scenario
      frontmatter, each citing a real platform-affordance dependency from
      MCP recon 2026-05-28 (not effort/complexity):
        - SR-01: function rows are not HTML5-draggable (Dart pointer-event
          DnD via _dndContext); the drag-drop leg uses the
          insertIntoCodeMirror END STATE via the CM6 view dispatch with
          usedFallback:true. The drag-drop flow is split to
          input-functions-ui.md (NOT in this scenario's owned
          ui_coverage_responsibility), so the JS-substituted end-state is
          an acceptable affordance-gap SR, not a pyramid-sub-rule violation.
        - SR-02: the canvas row->column mapping is non-linear and the
          literal Structure column is not among the visible canvas rows;
          the positive auto-bind type-match invariant is asserted FOR REAL
          against whichever Molecule/numeric column the probe selects (a
          name relaxation, not a console.warn downgrade).

      Open item for Gate F (not a Gate E failure): the chain YAML still
      lists this scenario's ui_coverage_responsibility as all THREE flows,
      whereas the scenario frontmatter lists only plus-icon with
      ui_coverage_split_to: [input-functions-ui.md]. Reconciling
      chain-level ownership vs the split is F-UI-COVERAGE-01's job; out of
      Gate E's scenario-paired scope.

      Proposal: accept the two existing scope_reductions[] entries (SR-01,
      SR-02) as-is and proceed. No spec change required at Gate E.
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-input-functions
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-powerpack-input-functions
    timestamp: 2026-05-20T00:00:00Z
    review_round: 1
    failure_keys: []
    claims:
      - check: A-STRUCT-MECH-01
        status: PASS
      - check: A-STRUCT-MECH-02
        status: PASS
      - check: A-STRUCT-MECH-03
        status: PASS
      - check: A-STRUCT-MECH-04
        status: PASS
      - check: A-STRUCT-MECH-05
        status: PASS
      - check: A-STRUCT-MECH-06
        status: PASS
      - check: A-STRUCT-03
        status: PASS
      - check: A-STRUCT-04
        status: PASS
      - check: A-LAYER-ALIGN-01
        status: PASS
      - check: A-CONT-01
        status: PASS
      - check: A-BUG-01
        status: PASS
      - check: A-MERIT-01
        status: PASS
      - check: A-MERIT-02
        status: PASS
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T12:33:30Z
    spec_runs:
      - spec: input-functions-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 61
        failure_keys: []
---

# Add New Column dialog — function insertion mechanisms (plus icon, drag-and-drop, auto-bound column parameter)

Medium-length UI scenario verifying the function-insertion mechanics of
the Add New Column dialog on the SPGI chem dataset: inserting a function
via the hover plus-icon, repeating via drag-and-drop, and verifying that
selected columns are auto-bound to the function's parameter when (and
only when) the column type matches the parameter type.

`pyramid_layer: ui-smoke` per `scenario-chains/powerpack.yaml` rev 1
(Rule 1 fallback — single source (SPGI), `classification: simple`, no
related bugs in frontmatter). NOT the chain's elected smoke witness;
`coverage_type: regression` is used here (the chain reserves
`coverage_type: smoke` for the elected smoke scenario at top-level
`add-new-column.md`).

`ui_coverage_delegated_to: add-new-column.md` — this scenario owns the
function-insertion-mechanism specialty flows (plus-icon, drag-n-drop,
auto-bound column parameter on type match) on its responsibility list
and delegates the basic dialog-open-and-add flow to the smoke witness.

`related_bugs: []` — no curated bug in `bug-library/powerpack.yaml`
directly intersects this scenario's reproduction surface. (Cross-cutting
candidates GROK-17109 and GROK-17004 affect the same `powerpack.dialogs.add-new-column`
+ `powerpack.dialogs.add-new-column-func` sub-features, but their
reproduction surfaces — datasync-persistence-and-recalc-on-rename and
complex-paste-handler-crash respectively — are not exercised by the
plus-icon / drag-and-drop / auto-bind flow here. See
`scenario-chains/powerpack.yaml :: bug_focused_candidates[]` — neither
GROK-17109 nor GROK-17004 lists this scenario in `spans`.)

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining (`depends_on: []` per chain rev 1).

## Scenarios

### Insert functions via + icon and drag-and-drop; auto-bind columns by type match

1. Open the `spgi.csv` dataset (the SPGI chem dataset that ships under
   `System:DemoFiles/chem/SPGI.csv`, containing a `Structure` column
   with chem-Molecule semType). **Verify:** the SPGI grid renders, with
   the `Structure` column showing molecule cell rendering.
2. Open the **Add New Column** dialog — click the `Add New Column`
   icon on the table view toolbar (or **Edit** > **Add New Column**
   from the top menu). **Verify:** the Add New Column dialog opens
   with its formula editor (CodeMirror), columns list (left), functions
   panel (right), and preview grid visible.
3. In the functions list on the right, hover the pointer over any
   function entry and click the `+` (plus) icon that appears on hover.
   **Verify:** the function is inserted into the formula text field
   with its parameter types as placeholders (e.g. hovering and clicking
   `+` on `Abs` inserts `Abs(num)` — function name plus its declared
   input parameter type token).
4. Clear the formula text field, then repeat adding the **same**
   function by dragging the function entry from the functions panel
   onto the formula text field. **Verify:** the same parameter-typed
   form appears in the formula text field (e.g. `Abs(num)`).
5. Click the `Structure` column in the columns list (left side). In
   the functions panel on the right, locate the `getCLogP` function,
   hover over it, and click the `+` icon. **Verify:** the function
   is inserted into the formula text field with the `Structure`
   column already passed as a parameter (e.g. `Chem:getCLogP(${Structure})`)
   — column auto-bound because its type (chem-Molecule) matches
   `getCLogP`'s expected input. **Verify:** the preview grid below
   reflects the function output computed against the `Structure`
   column.
6. Clear the formula text field, then repeat adding the **same**
   `getCLogP` function by dragging it from the functions panel onto
   the formula text field. **Verify:** the formula text field shows
   the same auto-bound form (`Chem:getCLogP(${Structure})`) and the
   preview grid again reflects the computed output.
7. Select any column of numeric type in the columns list. Add any
   function from the top of the functions list whose parameter is
   numeric (for example, `Abs`) — first by clicking its `+` icon,
   then by drag-and-drop. **Verify:** in both insertion modes, the
   selected numeric column is auto-bound as the function's parameter
   (e.g. `Abs(${SelectedNumericColumn})`), and the preview grid
   reflects the computed output.
8. Clear the formula text field. **Verify:** the formula text field
   is empty.
9. Click the `Id` column in the columns list. Click the sort-type
   icon in the top right corner of the dialog (the two-blue-arrows
   icon next to the functions panel header) and select sorting by
   `name`. **Verify:** the functions list is re-sorted alphabetically
   by function name.
10. Add the `Abs` function by drag-and-drop or by clicking its `+`
    icon. **Verify:** the function is inserted into the formula text
    field, but the `Id` column is **NOT** passed as a parameter —
    instead the parameter placeholder (`Abs(num)`) appears, because
    the `Id` column's type does not match `Abs`'s expected numeric
    input parameter type.

## Notes

- **target_layer: playwright** — chosen because a sibling
  `input-functions-spec.ts` already exists at the playwright layer
  (`existing-test-index.yaml` line 33501, layer `playwright`,
  category `PowerPack`). All other Add-New-Column siblings
  (`autocomplete-spec.ts`, `highlight-spec.ts`, `hints-spec.ts`,
  `functions-sorting-spec.ts`, `formula-refreshing-spec.ts`,
  `add-new-column-spec.ts`) are also Playwright specs, so the house
  style is consistent.
- **coverage_type: regression** — per the chain's smoke-scenario
  election (`ui_coverage_plan.smoke_scenario: add-new-column.md`),
  the top-level `add-new-column.md` owns the `smoke` slot. This
  narrower per-flow scenario lands in `regression`.
- **Helpers (already in registry, available for downstream
  Automator):** `softStep`, `loginToDatagrok`, `specTestOptions`,
  `stepErrors` from
  `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts` —
  used by the existing `input-functions-spec.ts`.
- **Bug-library status:** consulted —
  `bug-library/powerpack.yaml` exists and was scanned. No curated
  bug intersects this scenario's reproduction surface. Cross-cutting
  candidates GROK-17109 (calculated-column persistence across
  project save+datasync+reopen) and GROK-17004 (complex-paste-handler
  crash) touch the same `powerpack.dialogs.add-new-column` /
  `powerpack.dialogs.add-new-column-func` sub-features but their
  invariants are not exercised by the plus-icon / drag-and-drop /
  auto-bind flows. Chain `bug_focused_candidates[]` does not list
  this scenario in `spans` for either bug.
- **Decision log status:** queried — no
  `failed_attempts WHERE feature == powerpack` entries that touch
  the `AddNewColumn/input-functions.md` migration. No retry-skip
  approaches apply.
- **Atlas linkage (`derived_from:` provenance):**
  - `powerpack.dialogs.add-new-column` interaction "click function
    entry in functions panel → inserts function into formula with
    auto-bound matching column" is derived from
    `public/packages/PowerPack/src/tests/add-new-column.ts#L87`
    (atlas `interactions[]` entry, test-mined). This is the
    closest atlas-side analogue for the plus-icon / drag-and-drop
    / auto-bind mechanics this scenario exercises.
  - `powerpack.dialogs.add-new-column-func` interactions "click
    Add New Column toolbar icon → opens AddNewColumnDialog" and
    "Edit | Add New Column top-menu → opens AddNewColumnDialog"
    are code-derived from
    `public/packages/PowerPack/src/package.ts#L405`.
- **Source-text fixes applied:**
  - `add-new-column-dialog-title-case-normalization` — original
    body referred to the dialog as `*Add new column*` (sentence
    case in italics); migrated body uses the canonical product
    title `**Add New Column**` (title case in bold) to match the
    dialog's `[name="dialog-Add-New-Column"]` selector and the
    convention used in sibling migrated scenarios (`hints.md`,
    `highlight.md`).
  - `text-filed-to-text-field-typo` — original body contained the
    typo `text filed` (Steps 3 and 5); migrated body uses
    `text field`.
  - `fucntion-typo-corrected-to-function` — original body contained
    the typo `fucntion` (Step 10, `The fucntion should be added`);
    migrated body uses `function`.
  - `drag-n-drop-phrase-normalization` — original body alternated
    between `drag'n'drop` and `*plus* and drag'n'drop`; migrated
    body uses the normalized phrasing `drag-and-drop` and
    `+ icon` consistently. This matches the responsibility-key
    spellings `add-new-column-function-plus-icon` and
    `add-new-column-function-drag-drop` from the chain YAML.
- **Unresolved ambiguities (surfaced for downstream awareness):**
  - `sort-icon-glyph-described-not-named` — Step 9 describes the
    sort control as "the sort-type icon in the top right corner /
    two-blue-arrows icon". The existing `input-functions-spec.ts`
    selects this control via `[name="icon-sort-alt"]`. The original
    scenario text does not name the selector; the migrated body
    keeps the visual description for human readers and leaves
    selector-key resolution to the downstream Automator (the
    selector is already present in the on-disk spec, so no
    reference-file gap is introduced).
  - `id-column-type-classification` — Step 10 verifies that the
    `Id` column is NOT auto-bound to `Abs`'s numeric parameter,
    because the column type does not match. The on-disk
    `input-functions-spec.ts` line 244 asserts the formula text
    is `Abs(num)` (no substitution). The SPGI `Id` column is
    typed as string/qnum-like in the test fixture, so the
    type-mismatch claim holds; the migrated body preserves the
    original "type doesn't match the input parameter type"
    wording without further classification.
