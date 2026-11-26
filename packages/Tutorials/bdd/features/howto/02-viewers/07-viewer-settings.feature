@guide @tutorial
Feature: How-to · Viewers
  For the maintainer: the settings every viewer shares, on a scatter plot colored by RACE. The gear on the
  title bar opens the viewer's properties in the context panel, grouped in sections. Description: a
  title and a description typed into their rows show on the viewer. Legend: Legend Position moves the
  legend to the bottom. The Search box above the properties finds a setting in any section: "grid lines"
  leaves the two grid-line switches, and unchecking them clears the plot's background; "back color"
  leaves Style > Back Color, whose editor opens a palette (a light blue swatch). Nothing is saved
  as a default (the panel's own icons are not used); the scatter plot is added and demog opened off
  camera; nothing is made on the server. The Search box is named in howto-viewers.ts.

  Scenario: Change a viewer's settings
    The gear icon opens a viewer's settings: a title and description, the legend, and any setting found by name.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | Color | RACE |
    And scatter plot viewer takes 65% of the space it shares with the grid
    When the video says "Every viewer has its settings in the context panel"
    # caption: The gear on the title bar opens them
    And user hovers over empty plot space of scatter plot viewer
    And user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    # caption: Sections group the settings. Description: title and notes
    When user expands "Description" category in context panel
    And user enters "Height and weight" into "Title" property in context panel
    Then title of scatter plot viewer should have text "Height and weight"
    When user enters "5,850 patients of the demog study" into "Description" property in context panel
    Then description of scatter plot viewer should have text "5,850 patients of the demog study"
    # caption: Legend: where it goes
    When user expands "Legend" category in context panel
    And user selects "Bottom" in "Legend Position" property in context panel
    Then the legend of scatter plot viewer should be on the bottom
    # caption: Search finds a setting in any section
    When user clicks on context panel search
    And user types "grid lines" into context panel search
    And user unchecks "Show Vertical Grid Lines" property in context panel
    And user unchecks "Show Horizontal Grid Lines" property in context panel
    Then "Show Vertical Grid Lines" property of scatter plot viewer should be "false"
    And "Show Horizontal Grid Lines" property of scatter plot viewer should be "false"
    # caption: Style: colors, fonts and sizes
    When user clears context panel search
    And user types "back color" into context panel search
    And user clicks on editor of "Back Color" property in context panel
    And user picks the color "#9EDAE5" in the color picker
    And user clicks on empty plot space of scatter plot viewer
    Then the "view" area of scatter plot viewer should contain the color "#9EDAE5"
    And no errors should have been logged
