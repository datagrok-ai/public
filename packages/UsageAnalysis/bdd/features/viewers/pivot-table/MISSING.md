# Pivot table — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `pivot_viewer/CLAUDE.md` "Automation
surface" and `pivot-table.md`. The Table switch is in `pivot-table-table-switch.feature`.

## Needs a step

- **Scenario 6, the remembered aggregation type.** The `+` of the Aggregate row opens a column
  picker with an aggregation selector of its own (`pivot_grid.dart`, `showColumnComboBox`,
  `popupAggrSelector: true`); a pick there writes `PivotGrid.defaultAggrType`, which the
  `default aggregation` reading already reports. `user adds {string} to the {string} row of pivot
  table viewer` types the column only, and no step picks the type in that selector or cancels the
  picker without a pick. Needs `user adds {string} with {string} to the "aggregate" row of pivot
  table viewer` (the selector set before the column is committed, `bindings/tiers/viewers/pivot-table.ts`):

  ```gherkin
  Scenario: The + picker of the Aggregate row offers the type picked last, and Cancel adds nothing
    When user adds "HEIGHT" with "sum" to the "aggregate" row of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(AGE), sum(HEIGHT)"
    And the "default aggregation" reading of pivot table viewer should be "sum"
    When user clicks on the "remove aggregate chip sum(HEIGHT)" area of pivot table viewer
    And user clicks on the "add aggregate" area of pivot table viewer
    Then column picker popup should be visible
    When user clicks on status bar
    Then column picker popup should be absent
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    And the "default aggregation" reading of pivot table viewer should be "sum"
    And no errors should have been logged
  ```

  The reading is a static of `PivotGrid`, so the scenario must put it back to `avg` at its end.

## Needs a reading

- **The tag-editor header GROK-14995 duplicated** on a Table switch. `pivot-table-table-switch.feature`
  claims the Data row's chips; how many header rows the tag panel holds is not reported. Needs a
  `tag rows` reading (the panel's `grok-pivot-column-panel` count, `features/pivot_viewer_status.dart`)
  to claim `the "tag rows" reading of pivot table viewer should be 4` after the switch back.

## Product defect

- **Dragging a key column header inside the inner grid** is now a `@known-failure` scenario in
  `pivot-table-inner-grid-drag.feature` (core#467).

## Edge cases (now in `pivot-table-aggregation.feature`)

Written for the pivot table edge cases, now an `@edge-cases` scenario of its own, after the master scenarios of `pivot-table-aggregation.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Decided, steps requested**: a blank key clicked under Filtering Enabled is the `(empty)` category (an empty value is a category named `(empty)`), so the `filter label` should say `(empty)`. The claim needs a click on the blank key cell, which has no hit area yet (`key cell <row> of <column>`).

**Decided, steps requested**: a pivot on a column with thousands of categories (demog USUBJID) says that part of the categories was dropped (the silent-cut rule). Today `auto()` stops at 100 and the + picker has no cap, with no sign of either. The note has no wording yet, and the pivot reports no reading of it.

**Questions for owners**

None. ADD with nothing passing the filter adds an empty table, the user's choice (2026-10-30): claimed in `pivot-table-aggregation.feature`.

**Remark on the code docs**

- `pivot_viewer/CLAUDE.md` says a removed Group by column leaves `missingColumnError` over the grid. On the
  stand the removal re-picks a column (`look.auto` from `onColumnsRemoved`) and no error shows. The re-pick is
  the expected behaviour (a viewer that picks another column after its column is removed is correct), and it
  is requested from ApiTests below; the document is left as it is.

**API or Dart tests to write**

1. Dart regTest: the Data Manipulation crash test runs for `ViewerDescriptor.all.take(30)` and skips the pivot, the tile viewer and the points viewer.
2. Dart regTest: column names holding `, `, `[`, `]` or quotes through Save parameters and the history menu.
3. ApiTests: on demog-1000 grouped by DIS_POP, removing DIS_POP re-picks SEVERITY as Group by; the pivot shows no error, aggregates 5 rows and logs nothing.

Moved out of BDD (no user gesture); requested from ApiTests.

### Won't do

- A layout grouping by a column the table lacks: covered by the Dart regTest "Grouped by a column the table lacks" (`viewers_test.dart`).

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while folding the closed pivot table issues into `pivot-table-inner-grid.feature` and
`pivot-table-aggregation.feature`.

### Steps requested

- **The column picker of a tag row opens inside the window when the pivot is docked at the right edge
  (public#3690).** The picker the `+` of a tag row opens (`.d4-column-grid` with a
  `.d4-column-selector-backdrop`, `gestures.ts` `openColumnSelector`) has no noun, so no geometry claim
  can name it. Needs a `column picker popup` kind (the visible `.d4-column-grid:has(.d4-column-selector-backdrop)`),
  then:

  ```gherkin
  @realizes:public#3690
  Scenario: The Group by picker of a pivot docked at the right edge opens inside the window (public#3690)
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pivot table viewer
    When user docks pivot table viewer to the right edge of the view
    And user clicks on the "add group by" area of pivot table viewer
    Then column picker popup should be visible
    And column picker popup should be inside the window
    When user presses Escape
    Then column picker popup should be absent
    And no errors should have been logged
  ```
