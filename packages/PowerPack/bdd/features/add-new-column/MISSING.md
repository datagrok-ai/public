# Add New Column — what the features cannot say yet

Written during the calculated-columns gap round of 14 Oct 2026 (QA round of 13 Oct, help
`transform/add-new-column.md`) and reviewed on 18 Oct 2026 against the library's rules
(`public/libraries/bdd/CLAUDE.md`): a feature claims what only the browser shows, what an API or package
test can check goes there.

## Resolved

- **Rows that fail.** The gear next to the Type choice is named (`aria-label="Rows that fail"`), so
  `rows-that-fail.feature` opens it, sets Leave empty with an error column and Use value, and claims the
  preview's failure line, its Change... link, OK's balloon and one row each. Writing it found a bug, fixed
  in `add-new-column.ts`: the gear opened within the preview's debounce built a string editor for the value,
  and the type the preview then computed dropped the value on OK (package test "rows that fail: the value
  editor opened right after typing takes the type the formula returns").
- **core#806, core#784.** A double script input given whole numbers is fixed in core with its ddt test,
  so its known failure is gone from `script-in-formula.feature`; Append Tables now keeps OK disabled until a
  table is checked (a validator on the `tables` parameter, `tables_functions.dart`), which
  `append-tables.feature` claims instead of the NullError.
- **core#787, core#798.** Fixed in core with ddt tests (`calc_col_error_behavior_test.dart`): what a layout
  stores and what a function computes, so no scenario.
- **`formula error` scoped to the dialog** (`.add-new-column-dialog-root .cm-errort-div`): the Formula pane's
  line no longer satisfies it.
- **A script fixture named for a formula.** The stale-fixture sweep takes a `{time}` with no dash before it
  (`RUN_SUFFIX`, `src/runtime/server.ts`), so `BddDouble{time}` and `BddPowerSearch{time}` of a killed run go.
- **A known failure naming its logged error.** No longer needed: no known failure is left in this folder.
- **public#4293.** The preview ran the formula on a copy of the table still holding the edited column, so the
  result was "<name> (2)" on a dialog-made and an API-made column alike, once the 1 s preview had run (the first
  scenario cancels before it); `edit-column-formula.feature` claims it last, on a column made in the dialog.
- **public#4295.** String over a number adds the numbers as text, plain text the formula's text (the user's answer);
  `result-type-and-validation.feature` claims both, and an empty formula adding a plain empty column of the type.
- **An input whose nested input carries the label** (library, `byLabel`): a host input no longer goes by the label
  of an input nested in it, while a u2 `dart-input` bridge still goes by its inner root's (locator unit test).
- **A failed scenario leaving its dialog open over the next one.** Escape goes to the focused element, and the
  formula editor keeps it (by design: a stray Escape does not throw the typed formula away); the shell reset and a
  failed journey scenario now click the topmost dialog's close icon when Escape leaves it open (`closeOverlays`).

## Still open

Nothing planned.

## Won't do

- **A busy signal on the preview**: no claim reads "no error" right after a gesture; errors are claimed where
  they appear, and OK's result proves a valid formula.
- **An error column following a cell edit**: a recalculation, the package test's and ApiTests' "Rows the
  formula fails on".
- **A click right after OK reaching the context panel**: the freeze GROK-21024 keeps by design;
  `edit-column-formula.feature` makes the later columns current through the API.
- **What each formula computes on every row, the validator's cases, the error messages**: PowerPack package
  tests "Add new column" and ApiTests "DataFrame: Calculated columns".
- **`appendAsync`, `appendMergeAsync`, `grok.data.appendTables`, `recalculateFormulaColumns`; layout apply
  adding the formula columns in one batch (GROK-21114); UnixTimestamp/DateDiff precision (core#789,
  public#3851); data sync of `addNewCalculated` (core#788)**: what a function or a layout computes and stores.
- **Benchmarks** (calculated columns and layouts): the ApiTests benchmark category.
- **GROK-20931: the function at the caret as the current object, and the completion info beside the highlighted
  completion**: PowerPack package test "Add new column: function at caret becomes current object".
