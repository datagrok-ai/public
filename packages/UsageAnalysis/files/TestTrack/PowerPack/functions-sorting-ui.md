---
feature: powerpack
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  The dialog's columns list is a canvas-based grid with no addressable rows —
  a script cannot click a specific column, so these steps are verified by
  hand.
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


---
{
  "order": 5,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
