@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: a viewer can show fewer rows than its table without filtering the table. A scatter
  plot of demog-1000 (WEIGHT against HEIGHT: 872 rows have both) gets the Data > Filter formula
  ${AGE} > 44 and draws 466 rows while all 1,000 rows still pass the table's filter. Row Source set to
  Selected then limits it to the selected rows as well: nothing until rows are selected, and 8 for rows
  1 to 25 selected by their row headers in the grid (the 8 of them over 44). The formula is not cleared
  on camera: emptying the Filter row in the context panel leaves the plot filtered (reported). The plot
  is added and given most of the view off camera; demog-1000 is opened off camera; nothing is made on
  the server.

  Scenario: Show only some rows in a viewer
    A formula or the selection can limit the rows one viewer shows, while the table keeps them all.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    And scatter plot viewer takes 60% of the space it shares with the grid
    Then scatter plot viewer should show 872 rows
    When the video says "A viewer can show fewer rows than its table"
    # caption: The gear opens the viewer's settings
    And user hovers over empty plot space of scatter plot viewer
    And user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    # caption: Data > Filter: a formula picks the rows this viewer shows
    When user enters "${AGE} > 44" into "Filter" property in context panel
    Then scatter plot viewer should show 466 rows
    # caption: The grid keeps all 1,000 rows
    And all rows should pass the filter
    # caption: Row Source: here, only the selected rows
    When user selects "Selected" in "Row Source" property in context panel
    Then scatter plot viewer should show 0 rows
    # caption: Select rows in the grid: the plot shows those over 44
    When user selects rows 1 to 25 by their row headers
    Then 25 rows should be selected
    And scatter plot viewer should show 8 rows
    And no errors should have been logged
