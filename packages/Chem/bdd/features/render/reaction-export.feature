@journey @realizes:chem.reactions.export @eval-b7
Feature: The Export panel previews a reaction at print size and saves what it previews
  The Export panel, in place of two one-click buttons with the centre colours baked in and no choice of width, dpi or
  style, opens from the Reaction Sketcher's ribbon, a reaction cell's menu and the Reaction pane. It offers the presets
  ACS Document 1996 and Screen, the widths single column (8.3 cm), double column (17.8 cm) and free, and SVG, PNG at 300
  or 600 dpi, CDXML and Copy as picture. It offers no PDF. The centre colours are off for ACS 1996 and on for Screen.
  ACS 1996 in a single column is previewed 8.3 cm wide, and its PNG at 300 dpi is 980 px wide, black, with no legend.
  The two-step route of row 5 keeps print size there. The app's default route, whose agents are drawn as structures over
  its arrows, may be drawn smaller, and then the note says how much. What the panel draws for each choice (the rows, the
  legend, the SVG's size in points) is Chem's package tests ('Reactions: export').

  Background:
    Given user is logged in
    And the Chemical Engine is "Crux"

  Scenario: From the Reactions app's ribbon, on the default route
    Given user opens Chem's Reactions app on its default example
    And user watches downloads
    When user clicks Export... in the Reactions app's ribbon
    Then the Export panel should offer its presets, widths and formats, and no PDF
    And the Export panel's centre colours should be off for "ACS Document 1996" and on for "Screen"
    When user chooses "ACS Document 1996", "single column (8.3 cm)" and "PNG 300 dpi" in the Export panel
    Then the Export panel's preview should be 8.3 cm wide in the "publication" look
    And the Export panel's note should say how large the preview is against print size
    When user clicks Save in the Export panel
    Then a file "reaction.png" should have been downloaded
    And the saved PNG should be 980 px wide at 300 dpi, with no colour in it
    And no errors should have been logged

  Scenario: From a reaction cell's menu
    Given user opens test-reactions dataset
    And user watches downloads
    When user picks "Export reaction..." from the context menu of the "cell 5 of reaction" area of grid
    And user chooses "ACS Document 1996", "single column (8.3 cm)" and "PNG 300 dpi" in the Export panel
    Then the Export panel's preview should be 8.3 cm wide in the "publication" look
    And the Export panel's preview should be at print size
    When user clicks Save in the Export panel
    Then the saved PNG should be 980 px wide at 300 dpi, with no colour in it
    And no errors should have been logged

  Scenario: From the Reaction pane
    Given user opens test-reactions dataset
    And the context panel is open
    When user clicks on the "cell 5 of reaction" area of grid
    Then the context panel should show the current cell
    When user expands "Reaction" section
    And user clicks on "Export..." icon inside "Reaction" section
    And user chooses "ACS Document 1996", "double column (17.8 cm)" and "PNG 600 dpi" in the Export panel
    Then the Export panel's preview should be 17.8 cm wide in the "publication" look
    And the Export panel's preview should be at print size
    When user clicks Save in the Export panel
    Then the saved PNG should be 4205 px wide at 600 dpi, with no colour in it
    And no errors should have been logged
