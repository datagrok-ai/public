@viewers @realizes:GROK-17281
Feature: A layout saved before the legend position existed puts every legend on the right
  A layout saved on SPGI_v2 before viewers had a Legend Position (github #3203, GROK-17281) carries
  a histogram, a pie chart by Stereo Category, a scatter plot colored by it, a bar chart split by it
  and stacked by Primary Series Name, a line chart and a box plot, none of them with a
  legendPosition. Applied to SPGI, the five viewers that have a legend read Legend Position
  "Right". Four of them draw their legend on the right. The pie chart labels its slices and, with
  Legend Visibility Auto, draws no legend at this size; once its visibility is Always, it draws the
  legend on the right as well.
  The layout is opened from the computer through Open local file of the Browse toolbar: a file
  dropped on the view and a file opened that way both reach the platform's openFile, which applies a
  .layout file to the table view in front; the drop zone itself is not exercised. The fixture is the
  layout from the issue, unpacked into fixtures/nx.
  Not translated: the Context Panel as the place the position is read — the property is read from
  each viewer, which is what the panel shows.

  Scenario: Every viewer of the old layout has its legend on the right
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens spgi-3624 dataset
    And the browse panel is open
    When user uploads "fixtures/nx/spgi-legend-position-old.layout" through "Open local file" icon inside browse toolbar
    Then the open tableview should have 1 pie chart viewer
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 line chart viewer
    And the open tableview should have 1 box plot viewer
    And "legendPosition" property of pie chart viewer should be "Right"
    And "legendPosition" property of scatter plot viewer should be "Right"
    And "legendPosition" property of bar chart viewer should be "Right"
    And "legendPosition" property of line chart viewer should be "Right"
    And "legendPosition" property of box plot viewer should be "Right"
    And the legend of scatter plot viewer should be on the right
    And the legend of bar chart viewer should be on the right
    And the legend of line chart viewer should be on the right
    And the legend of box plot viewer should be on the right
    When user sets "legendVisibility" property of pie chart viewer to "Always"
    Then the legend of pie chart viewer should be on the right
    And no viewer of the current view should report an error
    And no errors should have been logged
    And no error or warning balloon should have been shown
