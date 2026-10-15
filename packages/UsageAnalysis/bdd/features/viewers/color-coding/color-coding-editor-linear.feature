@journey @viewers @realizes:viewers.grid
Feature: The color coding editor of a numeric column: type, Apply to, Min and Max, invert and Edit scheme
  The editor a column's header menu opens with Color Coding > Edit... ("Color-coding: AGE"), on
  demog-1000 (AGE: row 1 26, row 2 30, row 3 58). Each Color Coding item of the header menu replaces
  the column's coding, Off removes it (GROK-8282), and the items are a radio group with only the
  picked one checked (public#1455). In the editor, the Type choice switches the column between Off,
  Conditional and Linear; Apply to = text sets the column's text flag (that the renderer then paints
  the letters, and that background clears the flag, are in MISSING.md: no reading tells a text colour
  from a cell's, and no step claims the flag unset); Min and Max clamp the scheme, so two ages below Min take the same colour; the invert icon
  swaps the ends of the gradient; Edit scheme... opens "Edit color scheme", where Use absolute values
  (stops spread over the column's own range, so it differs only under a Min and Max) and removing a
  stop repaint the column at once (GROK-17736: the edited stops survive an invert).
  Translated from TestTrack Viewers/color-coding-ui.md steps 9-12 and the closed tickets public#1455,
  GROK-8282, GROK-17736. The claims are the colour the grid's renderer resolved for a cell, paired
  with the scheme or the type the column holds. Each scenario opens the editor itself and puts back
  what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then "AGE" column should have no color coding

  Scenario: Each Color Coding item of the header menu replaces the column's coding, and Off removes it
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    Then "AGE" column should be color-coded linearly
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user right-clicks on the "header AGE" area of grid
    And user hovers over "Color Coding" menu item
    Then "Color-Coding---Linear" menu item should be selected
    And "Color-Coding---Conditional" menu item should not be selected
    And "Color-Coding---Off" menu item should not be selected
    When user closes the context menu
    And user picks "Color Coding > Conditional" from the context menu of the "header AGE" area of grid
    Then "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user right-clicks on the "header AGE" area of grid
    And user hovers over "Color Coding" menu item
    Then "Color-Coding---Conditional" menu item should be selected
    And "Color-Coding---Linear" menu item should not be selected
    When user closes the context menu
    And user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should be the same
    And no errors should have been logged

  Scenario: The Type choice of the editor switches the column's coding
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    And Type input in "Color-coding: AGE" dialog should have value "Linear"
    When user selects "Conditional" in Type input in "Color-coding: AGE" dialog
    Then "AGE" column should be color-coded conditionally
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user selects "Off" in Type input in "Color-coding: AGE" dialog
    Then "AGE" column should have no color coding
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should be the same
    When user selects "Linear" in Type input in "Color-coding: AGE" dialog
    Then "AGE" column should be color-coded linearly
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: Apply to text sets the column's text flag
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user selects "text" in Apply to input in "Color-coding: AGE" dialog
    Then the text of "AGE" column should be color-coded
    When user selects "background" in Apply to input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: Min and Max clamp the scheme
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    And the "color of cell 1 of AGE" and "color of cell 2 of AGE" readings of grid should differ
    When user enters "50" into Min input in "Color-coding: AGE" dialog
    And user enters "60" into Max input in "Color-coding: AGE" dialog
    Then the "color of cell 1 of AGE" and "color of cell 2 of AGE" readings of grid should be the same
    And the "color of cell 1 of AGE" and "color of cell 3 of AGE" readings of grid should differ
    When user clears Min input in "Color-coding: AGE" dialog
    And user clears Max input in "Color-coding: AGE" dialog
    Then the "color of cell 1 of AGE" and "color of cell 2 of AGE" readings of grid should differ
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: The invert icon swaps the ends of the gradient
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    When user remembers the "color of cell 3 of AGE" reading of grid
    And user clicks on "exchange" icon in "Color-coding: AGE" dialog
    Then the color scheme of "AGE" column should be "#FF0000, #FFFFFF, #0000FF"
    And the "color of cell 3 of AGE" reading of grid should not be as remembered
    When user clicks on "exchange" icon in "Color-coding: AGE" dialog
    Then the color scheme of "AGE" column should be "#0000FF, #FFFFFF, #FF0000"
    And the "color of cell 3 of AGE" reading of grid should be as remembered
    When user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged

  Scenario: Edit scheme: absolute values ignore Min and Max, a stop is removed, and the edited scheme survives an invert
    When user picks "Color Coding > Linear" from the context menu of the "header AGE" area of grid
    And user picks "Color Coding > Edit..." from the context menu of the "header AGE" area of grid
    Then "Color-coding: AGE" dialog should be visible
    And the color scheme of "AGE" column should be "#0000FF, #FFFFFF, #FF0000"
    When user enters "50" into Min input in "Color-coding: AGE" dialog
    And user enters "60" into Max input in "Color-coding: AGE" dialog
    And user clicks on "Edit scheme..." action in "Color-coding: AGE" dialog
    Then "Edit color scheme" dialog should be visible
    When user remembers the "color of cell 3 of AGE" reading of grid
    And user checks "Use absolute values" input in "Edit color scheme" dialog
    Then "AGE" column should have tag ".%color-coding-linear-is-absolute" equal to "true"
    And the "color of cell 3 of AGE" reading of grid should not be as remembered
    When user unchecks "Use absolute values" input in "Edit color scheme" dialog
    Then the "color of cell 3 of AGE" reading of grid should be as remembered
    When user remembers the "color of cell 3 of AGE" reading of grid
    And user clicks on first Remove row button in "Edit color scheme" dialog
    Then the color scheme of "AGE" column should be "#FFFFFF, #FF0000"
    And the "color of cell 3 of AGE" reading of grid should not be as remembered
    When user clicks on OK button in "Edit color scheme" dialog
    Then "Edit color scheme" dialog should be hidden
    When user clicks on "exchange" icon in "Color-coding: AGE" dialog
    Then the color scheme of "AGE" column should be "#FF0000, #FFFFFF"
    When user clears Min input in "Color-coding: AGE" dialog
    And user clears Max input in "Color-coding: AGE" dialog
    And user clicks on CLOSE button in "Color-coding: AGE" dialog
    Then "Color-coding: AGE" dialog should be hidden
    And "AGE" column should be color-coded linearly
    When user picks "Color Coding > Off" from the context menu of the "header AGE" area of grid
    Then "AGE" column should have no color coding
    And no errors should have been logged
