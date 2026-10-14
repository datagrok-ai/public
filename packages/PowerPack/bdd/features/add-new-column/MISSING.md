# Add New Column — what the features cannot say yet

Written during the calculated-columns gap round of 14 Oct 2026 (QA round of 13 Oct, help
`transform/add-new-column.md`). Translation changed no code: every item below is either a step, a
name or a harness behaviour that is missing, with the scenario it would let us write, or a question
that needs a hand walk before it becomes a claim.

## 1. Rows that fail — the error behavior gear (blocked: the gear has no name)

The gear next to the Type choice (`initErrorBehavior` in `src/dialogs/add-new-column.ts`) opens the
"rows that fail" popup: **If a row fails** (Leave empty / Use value), **Value**, **Error column**,
**Column name**. It is an `ui.iconFA('cog')` with the class `add-new-column-error-behavior-icon` and
no `name` and no `aria-label`, so no phrase reaches it: the `icon` kind matches by aria, name or the
Dart `icon-{q}` name. The preview's failure line carries a **Change...** link that opens the same
popup, but that line is shown only once a setting is on — so the popup cannot be reached at all.

Needed, one of:

- a name on the gear in PowerPack (`name="icon-error-behavior"`, or an `aria-label` such as "Rows that
  fail"), which makes it `"error-behavior" icon in "Add New Column" dialog`;
- or a package element in `bindings/add-new-column.ts`:
  `element('error behavior icon', {selector: '.add-new-column-dialog-root .add-new-column-error-behavior-icon'})`.

And for the edit scenarios, a fixture step that makes a calculated column with an error behavior
through the API (the existing `user adds a calculated column {string} with formula {string}` passes
none), e.g. `user adds a calculated column {string} with formula {string} and an error column {string}`
→ `addNewCalculated(name, formula, {onError: {mode: 'empty', errorColumn}})`.

The popup's inputs are portaled outside the dialog; whether `"If a row fails" input` resolves without
a scope has not been probed (it could not be opened).

Scenarios, to be run and corrected once the gear has a name — the popup's phrases and the exact
counts and messages are taken from the source, not from a run (fixture: the three-row table of
`result-type-and-validation.feature` with a text column `d` = `2024-01-05`, `n/a`, `2024-03-07`):

```gherkin
  Scenario: A formula that fails on some rows: left empty, and the dialog says so
    When user clicks on "Add New Column..." icon
    And user types "parsed" into column name input
    And user types "DateParse(${d})" into formula editor
    And user clicks on "error-behavior" icon in "Add New Column" dialog
    And user selects "Leave empty" in "If a row fails" input
    And user checks "Error column" input
    Then "Column name" input should have value "parsed errors"
    And formula error should contain text "1 of 3 preview rows failed and will be left empty."
    And formula error should contain text "Messages are in \"parsed errors\""
    When user clicks on OK button in "Add New Column" dialog
    Then a warning balloon containing "1 of 3 rows failed and were left empty" should have been shown
    And the table should have a column "parsed errors"
    And the value of "parsed errors" column in row 2 should be "FormatException: 'n/a' is not a valid DateTime."
    And the value of "parsed" column in row 2 should be ""

  Scenario: Use value fills the rows that fail
    When user clicks on "Add New Column..." icon
    And user types "num" into column name input
    And user types "ParseInt(${d})" into formula editor
    And user selects "int" in column type input
    And user clicks on "error-behavior" icon in "Add New Column" dialog
    And user selects "Use value" in "If a row fails" input
    And user types "-1" into "Value" input
    Then formula error should contain text "filled with -1"
    When user clicks on OK button in "Add New Column" dialog
    Then the value of "num" column in row 2 should be "-1"

  Scenario: The error column follows a cell edit, and Change... reopens the setting
    When user double-clicks on the "cell 1 of d" area of grid
    And user presses Control+A in cell editor
    And user types "garbage" into cell editor
    And user presses Enter
    Then the value of "parsed errors" column in row 1 should be "FormatException: 'garbage' is not a valid DateTime."
    When user clicks on "Change..." link in "Add New Column" dialog
    Then "If a row fails" input should have value "Leave empty"

  Scenario: A renamed error column keeps its name when the layout is applied again (core#787)
    When user picks "Column Properties..." from the context menu of the "header parsed errors" area of grid
    And user types "why" into "New name:" input in "parsed errors" dialog
    And user clicks on OK button in "parsed errors" dialog
    And user saves the layout of the current table view
    And user loads the saved layout
    Then the table should have a column "why"
    And the table should not have a column "parsed errors"

  # core#798: auto type, Use value, every row fails -> a string column
  @known-failure
  Scenario: Use value on a formula that fails on every row keeps the computed type
    When user clicks on "Add New Column..." icon
    And user types "all bad" into column name input
    And user types "DateParse(${s})" into formula editor
    And user clicks on "error-behavior" icon in "Add New Column" dialog
    And user selects "Use value" in "If a row fails" input
    And user clicks on OK button in "Add New Column" dialog
    Then "all bad" column should have type "datetime"
```

The values of every row, the error messages and the behavior persisting in the column's tag are the
PowerPack package test "Add new column: rows that fail" and ApiTests "DataFrame: Calculated columns:
Rows the formula fails on" — the scenarios above claim the dialog and one row each.

## 2. public#4293 — the Edit Column Formula preview names the column "<name> (2)"

Ready as a known failure, but a known failure stops with the dialog open, and **Escape does not close
the dialog while the formula editor has the focus** (`codeMirrorDiv` keydown handler stops it, so the
editor's own popups can take it) — the dialog then stays over the next scenario of the journey, and
over the next feature at the shell reset. Needed: the harness closing a dialog a failed scenario left
open by its CANCEL (or close icon) when Escape leaves it on the page.

```gherkin
  # public#4293: the preview names the column under edit "<name> (2)"
  @known-failure:<line>
  Scenario: The preview of the dialog names the column by its own name
    Given the "next" column is the current object
    Then the context panel should show "next"
    Given Formula pane in context panel is expanded
    When user clicks on "Edit in dialog" button in context panel
    Then "Edit Column Formula" dialog should be visible
    And the preview grid should show "next" computed as numbers
    And the preview grid should not show "next (2)"
```

Not settled: in two runs here the claim passed — "next (2)" was not in the preview of a column
made through the API, while the bug reproduces by hand on a column made in the dialog (public#4293
screenshot). Either the column's origin matters or `should not show` read a preview before the
recalculation named the column; a probe of the preview's columns over time decides it before the
scenario goes in.

## 3. Library: a checkbox inside a multi-choice input

`"second" input in "Append Tables" dialog` resolves to two elements: the `second` bool input and the
whole **Tables** multi-choice host (`.ui-input-tables-multi-choice`), whose items carry the label
"second". `user checks` then sets the first checkbox of the first match — the host's `first` — so
`second` was never checked, and `should be checked` passed on the same host. The append ran over one
table. `append-tables.feature` uses the `checkbox` kind, which matches only the item. A fix in the
input kind: a host does not match by a label that belongs to an input nested inside it. Not filed as
an issue — for the library's owner to take or file.

## 3b. An existing feature that fails on the local stand (not touched)

`formula-refreshing.feature:127` (the Data sync round trip): after the project is reopened from
Browse > Dashboards, `user clicks on the "header Weight3" area of grid` leaves the reopened Project
current, and `the context panel should show "Weight3"` fails — 3 of 3 runs, alone too. With the
project opened through the API, a header click right after the open does reach the panel, so the
gallery open is the difference (the project seems to become current after the click). The feature is
not ours and is left as it is; `the "Weight3" column is the current object` before the pane claims
would make it pass, if the pick is not its subject.

## 4. Signals and steps

- **The dialog's validation has no "done" signal.** The type-mismatch and missing-column errors are
  written after the preview computes (1000 ms debounce, then `updateError`), so "no error under the
  editor" right after a gesture reads the empty line the dialog started with. The features claim an
  error only where it appears or goes away, and let OK's result prove a valid formula. Needed: a
  busy state on the preview while it recomputes (e.g. `aria-busy` from `setUpdateIndicator`), or a
  `formula error should settle` step that waits for it.
- **`formula error` is not scoped to the dialog.** The element is `.cm-errort-div`, which the Formula
  pane in the context panel uses too; an exact "" passes when any match is empty. The binding should
  be `.add-new-column-dialog-root .cm-errort-div` (the Formula pane would get its own element).
- **A known failure about a logged error cannot name the error.** `append-tables.feature` pins core#784
  to `no errors should have been logged`, so any other error on that path would count as the known
  one. Needed: `an error containing {string} should have been logged`, to pin it to
  `NullError: method not found: 'fold$2' on null`.
- **A script fixture named for a formula escapes the sweep of a killed run.** `script-in-formula.feature`
  names its script `BddDouble{time}`, since the formula calls it as an identifier; the stale-fixture
  sweep recognises only a `-<13 digits>` suffix (`RUN_SUFFIX`, `src/runtime/server.ts`), so a run killed
  before its feature end leaves the script behind (`power-search-results.feature`'s `BddPowerSearch{time}`
  has the same gap). Needed: a suffix form the sweep knows that is valid in an identifier, e.g. `_<ms>`.

## 5. Questions for a hand walk (not claimed)

1. **Plain text or string over a formula that computes a number is refused** — confirmed by hand,
   filed as [public#4295](https://github.com/datagrok-ai/public/issues/4295). `${x} + 1` with Type
   "plain text" (a text template: `Hello ${s}` gives "Hello a") or "string" shows "Result column
   type (int) doesn't match with current column type (string). Change column type or modify
   formula." and disables OK (`VALIDATION_TYPES_MAPPING` has no string entry), while a formula
   whose auto type is string (`${s} + ${s}`) is taken as a template ("a + a"). Intended? The legacy
   `help/transform/add-new-column-test.md` expects "Forcing int or string changes the preview".
2. **The context panel drops a click made right after the dialog's OK.** The OK makes the edited
   column current from code, and a cell clicked in the next second never reaches the panel (it stays
   on the edited column for good); under three workers a click right after a CANCEL was dropped as
   well. This is the freeze GROK-21024 kept by design. `edit-column-formula.feature` picks a column by
   its header once and then through `the "<col>" column is the current object`; a signal for the end
   of the freeze would let the later scenarios click as a person does.

## 6. Dropped by rule

- What each formula computes on every row, the validator's cases, the error messages: PowerPack
  package tests "Add new column" and ApiTests "DataFrame: Calculated columns".
- `appendAsync`, `appendMergeAsync`, `grok.data.appendTables`, `recalculateFormulaColumns`: ApiTests.
- Layout apply adding the formula columns in one batch and restoring grid width, visibility and order
  (GROK-21114), UnixTimestamp/DateDiff precision (core#789, public#3851), data sync of a column added
  with `addNewCalculated` (core#788, not planned): what a layout or a function stores and computes.
- Benchmarks (Calculated columns and layouts): ApiTests benchmark category.
