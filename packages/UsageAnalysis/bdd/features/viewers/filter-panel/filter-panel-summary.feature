@journey @viewers @realizes:viewers.filters
Feature: Filter panel indicator with viewers filtering
  The header counter of the panel counts the cards that restrict rows while a scatter plot zoom, a
  bar click, a pie slice, a trellis cell and a PC plot slider narrow the table further; its tooltip
  names each card with its categories or its range, and the panel's Reset filters icon empties the
  counter and gives back every row, the viewers' share included. One journey on demog-1000 (the
  md's spgi-100 with Stereo Category and Average Mass carried over as RACE and AGE, operator
  ruling D4): RACE Caucasian and AGE 30 to 60 leave 633 rows.
  Each viewer's click is shown to narrow the rows further than the viewers before it left them.
  The counter counts the panel's cards and nothing else — a bar click takes 633 rows to 303 and
  leaves it at 2, which is the behaviour the operator confirmed, and its tooltip names those cards
  with their categories and range. The md's "reflects all active filters" belongs elsewhere: the
  full list, the viewers' share included, is the tooltip of the "?" icon in the title bar of the
  Filters viewer — every source names itself there as "column: criterion", so a zoom reads
  "HEIGHT: [141.06,195.19]" and not "Scatter plot".
  Not translated: Scaffold Tree belongs to Chem (D3); AGE's range is set through the card's state.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a card for "RACE" to the filter panel
    And user clicks on the "category Caucasian of RACE" area of filter panel
    And user adds a range filter on "AGE" from 30 to 60
    Then 633 rows should pass the filter
    And counter of filter panel should have text "2"

  Scenario: Viewers filtering on top of the cards leave the counter at the cards
    When user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    And user drags a zoom box over the "view" area of scatter plot viewer
    Then fewer than 633 rows should pass the filter
    And counter of filter panel should be visible
    And counter of filter panel should have text "2"
    When user adds a bar chart viewer with:
      | Split | SEX |
    And user sets "On Click" property of bar chart viewer to "Filter"
    And user remembers the "rows shown" reading of filter panel
    And user clicks on the "bar M" area of bar chart viewer
    Then no rows where "SEX" is "F" should pass the filter
    And the "rows shown" reading of filter panel should be lower than remembered
    And counter of filter panel should have text "2"
    When user adds a pie chart viewer with:
      | Category | DIS_POP |
    And user sets "On Click" property of pie chart viewer to "Filter"
    And user remembers the "rows shown" reading of filter panel
    And user clicks on the "slice RA" area of pie chart viewer
    Then no rows where "DIS_POP" is "UC" should pass the filter
    And the "rows shown" reading of filter panel should be lower than remembered
    And the "rows shown" reading of filter panel should be at least 1
    And counter of filter panel should have text "2"
    When user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    And user remembers the "rows shown" reading of filter panel
    And user drags the max handle of the "AGE" range slider of pc plot viewer by 450 pixels
    Then the "filtering" reading of pc plot viewer should be "true"
    And the "rows shown" reading of filter panel should be lower than remembered
    And counter of filter panel should have text "2"
    When user hovers over counter of filter panel
    Then tooltip should contain the text "[30,60]"
    And tooltip should contain the text "RACE"
    And tooltip should contain the text "Caucasian"
    And tooltip should contain the text "AGE"
    And no errors should have been logged

  Scenario: A trellis cell filters on top of the cards and the counter stays at the cards
    When user adds a trellis plot viewer with:
      | X Column Names | SEVERITY     |
      | Y Column Names | SEX          |
      | Viewer Type    | Scatter plot |
    And user sets "On Click" property of trellis plot viewer to "Filter"
    And user remembers the "rows shown" reading of filter panel
    And user clicks on the "cell None | M" area of trellis plot viewer
    Then no rows where "SEVERITY" is "High" should pass the filter
    And the "rows shown" reading of filter panel should be lower than remembered
    And the "rows shown" reading of filter panel should be at least 1
    And counter of filter panel should have text "2"
    When user hovers over filter panel
    And user hovers over help icon of filter panel
    # the cell names both of the trellis axes in the summary, next to the cards
    Then tooltip should contain the text "SEVERITY: None"
    And tooltip should contain the text "RACE: Caucasian"
    And no errors should have been logged

  Scenario: The summary icon of the Filters title bar names the viewers' share as well as the cards
    When user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    And user drags a zoom box over the "view" area of scatter plot viewer
    And user adds a bar chart viewer with:
      | Split | SEX |
    And user sets "On Click" property of bar chart viewer to "Filter"
    And user clicks on the "bar M" area of bar chart viewer
    And user adds a pie chart viewer with:
      | Category | DIS_POP |
    And user sets "On Click" property of pie chart viewer to "Filter"
    And user clicks on the "slice RA" area of pie chart viewer
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    And user drags the max handle of the "AGE" range slider of pc plot viewer by 450 pixels
    And user hovers over filter panel
    Then help icon of filter panel should be visible
    When user hovers over help icon of filter panel
    # the two cards of the panel, each with the criterion it holds
    Then tooltip should contain the text "RACE: Caucasian"
    And tooltip should contain the text "AGE: [30,60]"
    # HEIGHT has no card: this line is the scatter plot's zoom and nothing else
    And tooltip should contain the text "HEIGHT: ["
    And tooltip should contain the text "SEX: M"
    And tooltip should contain the text "DIS_POP in [RA]"
    # the pc plot writes its ends with a space after the comma, where the card writes none
    And tooltip should contain the text "AGE: [18, "
    And tooltip should contain the text "Click for help (F1)"
    And no errors should have been logged

  Scenario: Reset filters empties the counter and gives back every row
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And counter of filter panel should be hidden
    And the "filters" reading of filter panel should be 0
    And the "cards" reading of filter panel should be "AGE, RACE"
    And the "filtering" reading of pc plot viewer should be "false"
    And no errors should have been logged
