@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: shaded ranges (formula annotation regions) on a histogram of AGE, demog-1000. The X axis's
  context menu, Annotations > Add Region, adds a region from the first to the third quartile (36 to 56: the middle
  half of the ages, 517 of the 1000 rows with the ties at both ends) and opens the Formula Lines dialog on it, where
  it is titled "Middle half". The same menu adds a second region, whose bounds are typed as formulas, ${AGE} = 65 and
  ${AGE} = 89 (the oldest patient is 89), titled "65 and over"; pointing at it above the bars (the histogram
  hit-tests a region only where no bin is under the pointer) shows its title and its 85 rows. Both stay on the
  viewer (its Annotation Regions property, in the browser only). The scatter plot's drawn regions are clip 02; this
  one is about a range of one column, which the histogram, bar chart, line chart and box plot axes all offer.
  demog-1000 is opened and the histogram added off camera; nothing is made on the server.

  Scenario: Shade a range of values on a chart
    Shade a range of an axis and name it: the middle half of the values, or any range you type.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a histogram viewer with:
      | Value | AGE |
    And histogram viewer takes 65% of the space it shares with the grid
    When the video says "Shade a range of values to point it out"
    # caption: Right-click the X axis: Annotations > Add Region
    And user picks "Annotations > Add Region" from the context menu of the "x axis" area of histogram viewer
    Then "Formula Lines" dialog should be visible
    # caption: It spans the first to the third quartile: the middle half
    When user enters "Middle half" into "Title" text input in "Formula Lines" dialog
    And user clicks on OK button in "Formula Lines" dialog
    Then histogram viewer should have a 'region "Middle half"' area
    # caption: Add another, then type its bounds: AGE from 65 to 89
    When user picks "Annotations > Add Region" from the context menu of the "x axis" area of histogram viewer
    Then "Formula Lines" dialog should be visible
    When user enters "${AGE} = 65" into "Formula 1" text area in "Formula Lines" dialog
    And user enters "${AGE} = 89" into "Formula 2" text area in "Formula Lines" dialog
    And user enters "65 and over" into "Title" text input in "Formula Lines" dialog
    And user clicks on OK button in "Formula Lines" dialog
    Then histogram viewer should have a 'region "65 and over"' area
    And histogram viewer should have a 'region "Middle half"' area
    # caption: Point at a region: its name and how many rows fall in it
    When user hovers over the 'region "65 and over"' area of histogram viewer
    Then tooltip should contain text "85 rows"
    And no errors should have been logged
