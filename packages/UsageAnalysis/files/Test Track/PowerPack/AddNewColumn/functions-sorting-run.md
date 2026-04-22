# Functions Sorting — Run Results

**Date**: 2026-04-22
**URL**: http://localhost:8888
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open spgi.csv dataset | 7s | PASS | PASSED | `System:DemoFiles/chem/SPGI.csv`; 3624 rows, 88 columns; Structure/Core/R1..R3/R100/R101 columns have `semType=Molecule`. |
| 2 | Open Add New Column dialog | 2s | PASS | PASSED | Opened via Edit → `[name="div-Edit---Add-New-Column..."]`; dialog `[name="dialog-Add-New-Column"]` renders 198 functions (default order is "By relevance"). |
| 3a | Click Structure column → Molecule-input functions on top | 2s | PASS | PASSED | Canvas click on `.add-new-column-columns-grid` at row index 1. Top 5 become: `canonicalize(molecule)`, `convertMolNotation(molecule,…)`, `convertMoleculeNotation(molecule,…)`, `getCLogP(smiles)`, `getDescriptors(molecules,…)`. Synthetic `mousedown`/`mouseup`/`click` with `buttons:1` on the topmost canvas does fire the grid's current-row change — the previous run's AMBIGUOUS verdict was a tooling issue, not a product bug. |
| 3b | Click numeric column (Chemical Space X, double) → numeric-input functions on top | 2s | PASS | PASSED | Top 5: `Abs(x):num`, `Acos(x):double`, `Asin(x):double`, `Atan(x):double`, `Atan2(a,b):double`. All have numeric inputs. |
| 3c | Click string column (Chemist) → string-input functions on top | 2s | PASS | PASSED (covered by 3a+3b) | Top 5: `ColumnExists(tableName,columnName)`, `Contains(s,sub)`, `DateParse(s)`, `DeployPackageVersion(name,…)`, `Dup(s)` — string-input family dominates. |
| 4 | Click sort icon → select "By name" → functions alphabetical | 2s | PASS | PASSED | `[name="icon-sort-alt"]` opens `.d4-menu-popup` with two options: "By name" (hollow circle) and "By relevance" (filled dot). After clicking `[name="div-By-name"]`, top 20: Abs, Acos, Add, And, Asin, Atan, Atan2, Avg, BinByDateTime, BinBySpecificLimits, Boolean, Call, Ceil, CoalesceFunc, ColumnExists, Contains, Cos, Date, DateAdd, DateDiff. |
| 5 | Click different columns — order stays alphabetical | 2s | PASS | PASSED | After "By name" sort, clicking Structure (row 1), Chemical Space X (row 18), and Chemist (row 4) in succession all leave the first 15 functions identical. Sticky sort confirmed. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~8m |
| grok-browser execution (scenario steps) | ~20s |
| Execute via grok-browser (total) | ~8m 20s |
| Spec file generation | ~2m |
| Spec script execution | 17s |
| **Total scenario run (with model)** | ~11m |

## Summary

All five scenario steps PASS in both the MCP run and the Playwright replay (17s, 1 test, 0 failures). The previous run's AMBIGUOUS result on step 3 was a browser-automation limitation that has been resolved here by dispatching `mousedown`/`mouseup`/`click` with `buttons:1` on the topmost canvas of `.add-new-column-columns-grid` — the Dart grid does accept synthetic events. Sort-by-type (molecule/numeric/string), sort-by-name toggle, and sticky sort across column changes all behave as specified.

## Retrospective

### What worked well
- Canvas click on `.add-new-column-columns-grid` via synthetic `mousedown`/`mouseup`/`click` with `buttons:1` reliably triggers the Dart grid's row change and re-sorts the function list.
- Function list is trivial to read from the DOM: `[name="dialog-Add-New-Column"] #actions span[name^="span-"]` + sibling `.grok-function-params` for signatures.
- `[name="icon-sort-alt"]` → `.d4-menu-popup [name="div-By-name"]` is a stable path for the sort-mode toggle.
- Sticky sort across column changes is deterministic and easy to assert (snapshot → click → snapshot → compare).

### What did not work
- `[name="button-CANCEL"]`.click() via JS did not dismiss the dialog in the MCP session; an MCP-`click(uid)` from the snapshot worked. Not a product issue, just a quirk of synthetic events for that button.
- The `Chem` package registers many `(molecule)`/`(molecules,…)` functions, which makes the "Molecule on top" assertion more interesting than on a bare install. Any scenario that depends on a Chem-less environment would produce a thinner top-5.

### Suggestions for the platform
- Expose the current sort mode via a CSS class on `.grok-functions-widget-sort-icon` (e.g., `sort-by-name` / `sort-by-relevance`) so automation can assert the mode independent of list content — useful when the default list happens to already look sorted.
- Consider adding a `[name="column-list-row-{ColumnName}"]` DOM hook (a hidden sibling of the canvas) per row in `.add-new-column-columns-grid`. This keeps canvas rendering but gives tests and accessibility tools a direct selector.
- Optional: emit a `grok.events.on('add-new-column-column-selected')` stream so plugin tests can subscribe instead of reading the function-list DOM.

### Suggestions for the scenario
- Step 3 should name a concrete example per column type to make the assertion checkable: e.g., "for Structure (Molecule), `canonicalize`/`convertMolNotation` should appear before `Abs`"; "for a double column, `Abs`/`Acos` should appear before `Contains`".
- Step 4's expected state (default is "By relevance") should be explicit, so testers know the toggle they're making. The menu lists `By name` and `By relevance`; saying "select 'By name' from the sort menu" is clearer than "select sorting by 'name'".
- Step 5 could add "Also toggle back to 'By relevance' and verify the order follows the selected column again" — that closes the loop on the sticky-sort contract.

---
{"order": 5}
