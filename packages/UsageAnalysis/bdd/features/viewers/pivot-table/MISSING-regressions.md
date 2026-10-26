# Pivot table — closed issues the regression round of 26 Oct 2026 could not claim

Kept apart from `MISSING.md` (the regression round of closed issues; the edge-case round keeps its own notes there).

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed pivot table issues into `pivot-table-regressions.feature`.

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
