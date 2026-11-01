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
