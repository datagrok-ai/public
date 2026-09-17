@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot row source
  Which rows the cells are drawn from. Row Source and On Click cannot both be steering the filter,
  so the platform corrects whichever was set last: On Click = Filter moves Row Source to All, and
  Row Source = Filtered moves On Click back to None. Then the whole ladder of eight row sources,
  each read as the number of rows the trellis has and the number of cells that painted nothing —
  with a filter card on SEX = F and every Caucasian row selected, no two rungs read alike.

  Packing is off throughout, so the grid stays 2 by 4 whatever the row source leaves: a rung is
  read on the rows, not on a grid that reshaped itself.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names  | SEX          |
      | Y Column Names  | RACE         |
      | Viewer Type     | Scatter plot |
      | Pack Categories | false        |
    Then the "cells" reading of trellis plot viewer should be 8
    And trellis plot viewer should show 1000 rows

  Scenario: On Click Filter moves Row Source to All, and Filtered moves On Click back to None
    Then "Row Source" property of trellis plot viewer should be "Filtered"
    And "On Click" property of trellis plot viewer should be "None"
    When user sets "On Click" property of trellis plot viewer to "Filter"
    Then "Row Source" property of trellis plot viewer should be "All"
    And "On Click" property of trellis plot viewer should be "Filter"
    When user sets "Row Source" property of trellis plot viewer to "Filtered"
    Then "On Click" property of trellis plot viewer should be "None"
    And "Row Source" property of trellis plot viewer should be "Filtered"
    And no errors should have been logged

  Scenario: The filter and the selection every rung is read against
    When user adds a categorical filter on "SEX" keeping "F"
    Then 553 rows should pass the filter
    When user selects rows where "RACE" is "Caucasian"
    Then 896 rows should be selected
    When user sets "Row Source" property of trellis plot viewer to "All"
    Then trellis plot viewer should show 1000 rows
    And the "cells" reading of trellis plot viewer should be 8
    And the "blank cells" reading of trellis plot viewer should be 0
    And no errors should have been logged

  Scenario Outline: Row Source <source> feeds the cells <shown> rows
    When user sets "Row Source" property of trellis plot viewer to "<source>"
    Then "Row Source" property of trellis plot viewer should be "<source>"
    And trellis plot viewer should show <shown> rows
    And the cells of trellis plot viewer should be 2 wide and 4 tall
    And the "cells drawn" reading of trellis plot viewer should be 8
    And the "blank cells" reading of trellis plot viewer should be <blank>
    And the "distinct cell signatures" reading of trellis plot viewer should be <pictures>
    And "On Click" property of trellis plot viewer should be "None"
    And no errors should have been logged

    Examples:
      | source            | shown | blank | pictures |
      | All               | 1000  | 0     | 8        |
      | Filtered          | 553   | 4     | 5        |
      | Selected          | 896   | 6     | 3        |
      | FilteredSelected  | 480   | 7     | 2        |
      | SelectedOrCurrent | 896   | 6     | 3        |
      | CurrentRow        | 1     | 7     | 2        |
      | MouseOverGroup    | 0     | 8     | 1        |
      | MouseOverRow      | 0     | 8     | 1        |

  Scenario: Putting the ladder back
    When user sets "Row Source" property of trellis plot viewer to "Filtered"
    And user clears the row selection
    Then no rows should be selected
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then "SEX" filter card should be absent
    And all rows should pass the filter
    And trellis plot viewer should show 1000 rows
    And the "blank cells" reading of trellis plot viewer should be 0
    And no errors should have been logged
