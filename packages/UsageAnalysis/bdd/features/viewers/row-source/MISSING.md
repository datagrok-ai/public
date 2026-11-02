# Row Source edge cases — what the features cannot say yet

Written for the Row Source edge cases (now in `row-source-mouse-over-group.feature`): the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

## Still open

**Decided**

1. A viewer Filter `${AGE} > 44` follows a rename of AGE: core#1213, still open, requested from ApiTests in `../MISSING.md` (today the formula keeps the old name, is treated as empty and the viewer widens to every row, `validateColumnFormulaFields`).
2. Under MouseOverGroup, closing the viewer whose group is hovered (a pie) clears the group (a filter left by a closed viewer goes with it): written in `row-source-mouse-over-group.feature` on the `rows shown` of a MouseOverGroup scatter plot, core#1227, a known failure (the scatter plot keeps drawing the closed pie's 5 Asian rows its own `${AGE} > 44` keeps).

**Questions for owners**

None. A viewer Filter whose column is removed is dropped without a message, as the rules say (2026-10-30).

**API or Dart tests to write**

1. Dart regTest: two viewers sharing one Filter formula, one closed while the formula column is computing (`dataFrame.temp` flag).
2. Dart regTest: the Data Manipulation crash test looped over the row sources (it runs on Filtered only).
3. ApiTests: a Filter `${AGE} + 1` on a scatter plot (X = AGE), a histogram (Value = AGE) and a pie chart (Category = AGE) shows the "Invalid formula" error balloon and leaves 1000 rows shown; `${AGE} > 44` then shows 519.
4. ApiTests: a scatter plot under Row Source FilteredSelected with the 447 M rows selected and only F passing the filter shows 0 rows and logs nothing.

Moved out of BDD (no user gesture); requested from ApiTests.

## Won't do

- Row Source and Filter through a saved layout for the seven viewers: what a layout stores is an API test, and the viewers' own persistence features restore their settings.
