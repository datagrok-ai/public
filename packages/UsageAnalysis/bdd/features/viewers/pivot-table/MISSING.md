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

## Candidate product defect, not a gap

- **Dragging a key column header inside the inner grid** (`pivot-table-ui.md`, re-sort the rows):
  `user drags the "grid header DIS_POP" area of pivot table viewer to the "grid header None avg(AGE)"
  area` logs `Cannot read properties of null (reading '$lt')` (measured 23 Oct; the inner-grid feature
  describes the NullError from `col2screen` on every repaint). No ticket found in `datagrok-ai/core`
  or `public`; not a `@known-failure` until one exists.
