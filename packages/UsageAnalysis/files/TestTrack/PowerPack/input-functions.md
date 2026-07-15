---
feature: powerpack
target_layer: playwright
coverage_type: smoke
priority: p1
realizes_atlas: []
realizes: []
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

UI scenario verifying the function-insertion mechanics of the Add New
Column dialog on the SPGI chem dataset: inserting a function via the
hover plus-icon, repeating via drag-and-drop, and verifying that
selected columns are auto-bound to the function's parameter when (and
only when) the column type matches the parameter type.

The dialog's basic open-and-add flow is covered by `add-new-column.md`.

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining is required.

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

