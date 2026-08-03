---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
target_layer: manual-only
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-function-drag-drop
  - add-new-column-auto-bound-column-parameter
ui_coverage_delegated_to: add-new-column.md
companion_to: input-functions.md
produced_from: ui-affordance-split
split_date: 2026-05-27
split_rationale: |
  Two owned ui-smoke flows on the parent `input-functions.md` proved
  empirically unscriptable from outside the Dart runtime per MCP recon
  (cycle 2026-05-26-powerpack-automate-03, Round-2): drag-drop of a
  function entry onto the formula editor, and selection of a column row
  in the dialog's column-grid widget (private `columnsDf`). All
  JS-accessible paths refuted (DG.Grid.fromRoot stub, columnsDf not in
  grok.shell.tables, DG.Dialog wrapper carries no columnsDf, source-df
  currentCol does not propagate to dialog selectedColumn, synthetic
  mouse events absorbed by overlay canvas, keyboard arrow events
  Dart-event-bus-bound). The parent spec had to fall back to JS-API
  substitution (`document.execCommand('insertText', ...)` + source-df
  `currentCol` setter with `usedFallback: true` annotation), which
  violates the pyramid_layer: ui-smoke rule "JS API substitution NOT
  permitted for flows in ui_coverage_responsibility". This companion
  carries those affordances as manual-only coverage; the parent retains
  the `add-new-column-function-plus-icon` flow (trusted `.click()` on
  `[name="icon-plus"]` — Gate B PASSED).
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
---

## Setup

A clean Datagrok session. Open the `spgi.csv` dataset (the SPGI chem
dataset that ships under `System:DemoFiles/chem/SPGI.csv`, containing
a `Structure` column with chem-Molecule semType). Open the **Add New
Column** dialog — click the `Add New Column` icon on the table view
toolbar, or **Edit** > **Add New Column** from the top menu.

## Scenarios

### Drag-and-drop a function entry onto the formula editor

Verifies the drag-and-drop affordance for inserting a function from the
functions panel into the formula text field. The drop end-state must
match the plus-icon insertion end-state (parameter-typed form, e.g.
`Abs(num)`).

1. In the functions list on the right, locate any function entry (for
   example, `Abs`).
2. Drag the function entry from the functions panel onto the formula
   text field of the dialog editor. **Verify:** the formula text field
   shows the parameter-typed form (e.g. `Abs(num)`) — function name
   plus its declared input parameter type token.

### Auto-bound column parameter on type match (column selected in column-grid)

Verifies that selecting a column in the dialog's column-grid widget
auto-binds it as the function's parameter when the column type matches
the parameter type. The trigger is column-grid selection; the
behavioural assertion is the `${ColName}` form in the editor.

1. In the columns list on the left of the dialog, click the `Structure`
   column row.
2. In the functions panel on the right, locate the `getCLogP` function,
   hover over it, and click the `+` (plus) icon that appears on hover.
   **Verify:** the function is inserted into the formula text field
   with the `Structure` column already passed as a parameter
   (`Chem:getCLogP(${Structure})`) — column auto-bound because its type
   (chem-Molecule) matches `getCLogP`'s expected input.
3. Verify the preview grid below the editor reflects the function
   output computed against the `Structure` column.
4. Clear the formula text field, then repeat adding `getCLogP` by
   dragging it from the functions panel onto the formula text field.
   **Verify:** the formula text field shows the same auto-bound form
   (`Chem:getCLogP(${Structure})`) and the preview grid again reflects
   the computed output.
5. Select any column of numeric type in the columns list. Add any
   function whose parameter is numeric (for example, `Abs`) — first by
   clicking its `+` icon, then by drag-and-drop. **Verify:** in both
   insertion modes, the selected numeric column is auto-bound as the
   function's parameter (e.g. `Abs(${SelectedNumericColumn})`), and
   the preview grid reflects the computed output.

### Auto-bound NEGATIVE contract on type mismatch

Verifies that selecting a column whose type does NOT match the
function's parameter type results in the parameter placeholder
(parameter-typed form), not an auto-bound `${ColName}`.

1. Clear the formula text field.
2. Click the `Id` column in the columns list. Add the `Abs` function by
   drag-and-drop or by clicking its `+` icon. **Verify:** the function
   is inserted into the formula text field, but the `Id` column is
   **NOT** passed as a parameter — instead the parameter placeholder
   (`Abs(num)`) appears, because the `Id` column's type does not match
   `Abs`'s expected numeric input parameter type.

## Notes

- **Why manual-only**: MCP recon 2026-05-26 (cycle
  2026-05-26-powerpack-automate-03, Round-2) conclusively refuted all
  JS-accessible paths to the dialog's private column-grid `columnsDf`
  and the CodeMirror `cmView` reference. Playwright's `dragTo` on the
  function-entry → formula-editor pair did not trigger the PowerPack
  `ui.makeDroppable` subscription (seeded only by source-side Dart
  `dragstart` handler). Both affordances are exercised correctly by a
  human tester but not from outside the Dart runtime under current
  build. See parent scenario's `gate_verdicts.e.scope_reduction_proposal`
  for the full MCP recon evidence.
- **Parent scenario coverage**: `input-functions.md` retains ownership
  of `add-new-column-function-plus-icon` (the plus-icon insertion path,
  which uses trusted `.click()` on `[name="icon-plus"]` and is fully
  scriptable; Gate B PASSED on this in
  cycle-1 of 2026-05-26-powerpack-automate-03).
- **Reinstate as Playwright**: if a future PowerPack build exposes
  either (a) the column-grid `columnsDf` via `DG.Dialog.getOpenDialogs()`
  with a public accessor, OR (b) the CodeMirror `cmView` on the
  `.cm-editor` host, OR (c) a public `setSelectedColumn(name)` API on
  the dialog wrapper, the affordances become scriptable and the flows
  can be lifted back to the parent Playwright spec. The split is
  reversible; track via `companion_to` ↔ `ui_coverage_split_to`
  back-links.
- **Bug-library status**: the affordance gaps are not catalogued in
  `bug-library/powerpack.yaml` because they are platform-architecture
  constraints (Dart event-bus binding, private wrapper internals), not
  fixable bugs. If a JIRA ticket is filed to expose the affordances
  publicly, link it via `related_bugs[]` here.

---
{
  "order": 99
}
