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

## Edge cases (`pivot-table-edge-cases.feature`)

Written for `pivot-table-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Questions for owners**

1. Row Source Filtered with nothing passing, then ADD: hide ADD, or publish an empty table?
2. Group by a column with blanks and a click on the blank key cell under Filtering Enabled: what does the `filter label` say for the blank value?
3. Pivot on a column with thousands of categories (demog SUBJ): `auto()` stops at 100, the + picker does not. A cap or a warning?
4. `pivot_viewer/CLAUDE.md` says a removed Group by column leaves `missingColumnError` over the grid; on the stand the removal re-picks the column (`look.auto` from `onColumnsRemoved`) and no error shows. The feature claims the re-pick; the document should say which applies to a removal and which to a layout from another table.

**API or Dart tests to write**

1. Dart regTest: the Data Manipulation crash test runs for `ViewerDescriptor.all.take(30)` and skips the pivot, the tile viewer and the points viewer.
2. Dart regTest: column names holding `, `, `[`, `]` or quotes through Save parameters and the history menu.

### Won't do

- A layout grouping by a column the table lacks: covered by the Dart regTest "Grouped by a column the table lacks" (`viewers_test.dart`).
