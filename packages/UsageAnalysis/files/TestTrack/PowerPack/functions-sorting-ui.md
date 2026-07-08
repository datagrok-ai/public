---
feature: powerpack
target_layer: manual-only
pyramid_layer: ui-smoke
ui_companion_of: functions-sorting.md
companion_steps: [3, 4]
ui_coverage_responsibility:
  - add-new-column-functions-sort-by-type
manual_execution_notes: |
  The functions-panel sort-BY-TYPE flow (clicking a column in the Add New
  Column dialog's columns list so functions whose first input parameter
  matches the selected column's type rise to the top) is NOT automatable
  from a Playwright spec. The columns list is a canvas-based popup-mode
  DG.ColumnGrid (PowerPack/src/dialogs/add-new-column.ts:1079,
  widgetMode: true) with ZERO DOM-addressable per-row elements, and the
  Dart-side hit-test pipeline for this widget does NOT respond to
  synthetic DOM pointer/mouse events dispatched from a Playwright spec.
  Live MCP recon (dev.datagrok.ai, 2026-05-28) confirmed across four
  trigger approaches — synthetic MouseEvent triple-sequence, full
  pointer+mouse event chain, JS-API DG.Grid.fromRoot().dataFrame (null),
  and search-input filter + keyboard selection — that NONE cause the
  functions list to re-sort. A real user's OS-level mouse click DOES
  trigger the re-sort (the feature works for humans), but that path is
  unreachable from headless automation. These two steps are therefore
  carved out for manual visual verification; the canonical
  functions-sorting.md / functions-sorting-spec.ts retains the
  automatable flows (dialog open, sort-by-name, alphabetical ordering).
---

# Add New Column dialog — functions-panel sort-BY-TYPE (manual)

Manual companion to `functions-sorting.md`. Covers the two scenario steps
whose trigger (selecting a column in the dialog's canvas-based columns
list) cannot be automated. Verify these by hand with a real mouse on a
live Datagrok instance.

## Preconditions

1. Open the SPGI dataset (`System:DemoFiles/chem/SPGI.csv`). The grid
   renders with its chem `Structure` column (semType `Molecule`) plus
   numeric and string columns.
2. Open the **Add New Column** dialog — click the `Add New Column` icon on
   the table view toolbar (or **Edit** > **Add New Column**). The dialog
   opens with its formula editor, columns list, functions panel, and
   preview grid.
3. Ensure the functions panel sort mode is **"By relevance"** (click the
   sort icon `[name="icon-sort-alt"]` at the top-right of the functions
   panel and select "By relevance" if it is currently "By name").

## Manual steps

### Step 3 (manual) — Structure column → Molecule-input functions on top

With a real mouse, click the `Structure` column row in the dialog's
columns list. **Verify:** the functions list on the right re-sorts so that
functions with `Molecule`-type input parameters appear at the top (for
example, `canonicalize(molecule)`, `convertMolNotation(molecule, ...)`,
`getCLogP(smiles)`, `getDescriptors(molecules, ...)`) ahead of functions
whose first input is numeric or string.

Note: the exact top-of-list function names depend on the server's current
function catalogue and may differ from the examples above. The contract
being verified is that the functions list **re-orders** on the column
click and that a Molecule-input family is brought toward the top.

### Step 4 (manual) — numeric column → numeric-input functions on top

With a real mouse, click a numeric column (for example `Chemical Space X`,
type `double`) in the columns list. **Verify:** the functions list
re-sorts so functions whose first input parameter is numeric appear at the
top (for example, `Abs(x)`, `Acos(x)`, `Asin(x)`, `Atan(x)`,
`Atan2(a, b)`) ahead of chem-input or string-input families. Repeat for
any other column type as available and verify the matching-parameter
family is on top in each case.

## Notes

- **Why manual.** The dialog's columns list is a canvas-based grid
  with no addressable rows, so a script cannot reliably click a
  specific column in it; a real user's mouse click works fine. This
  step is therefore verified by hand.
- **Automated sibling.** `functions-sorting.md` covers the dialog-open
  flow and the sort-icon → "By name" → alphabetical-ordering flow,
  both of which are reliably scriptable. The sticky-sort check (once
  "By name" is active, clicking columns must not change the order) is
  verified manually here as a companion check to the two manual steps
  above.
