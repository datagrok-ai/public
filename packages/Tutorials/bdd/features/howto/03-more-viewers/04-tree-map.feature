@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Two views of one hierarchy. The tree map from the toolbox
  (the toolbox then hidden with its sidebar icon, and the map dragged by its title bar to the right edge, on
  camera, so its selectors are not at the window's edge): SEX as the top level and
  DIS_POP as the second, WEIGHT as the color, all with the selectors on the map, and a rectangle clicked (its 1941
  rows selected). Then the sunburst from the Add viewer gallery (Charts package), its Hierarchy set to the same SEX
  and DIS_POP in the column dialog of its settings (the dialog keeps the order of its rows), and a segment clicked
  (820 rows). demog is opened off camera; nothing is made on the server.

  Scenario: Show a hierarchy with a tree map or a sunburst
    Categories inside categories: a tree map nests rectangles, a sunburst nests rings.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A tree map shows categories inside categories"
    And user clicks on tree map icon on toolbox
    Then tree map viewer should be added to the open tableview
    # caption: Hide the toolbox to give the map room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Drag its title bar to dock it on the right
    When user docks tree map viewer to the right edge of the view
    Then tree map viewer should be docked along the right edge of the view
    # caption: Pick the top level here
    When user selects "SEX" in first tree map split selector
    # caption: and add the next level beside it
    And user selects "DIS_POP" in second tree map split selector
    Then the "levels" reading of tree map viewer should be 2
    # caption: Color the rectangles by the average of a column
    When user picks "WEIGHT" in the "color" column selector of tree map viewer
    Then the "color column" reading of tree map viewer should be "WEIGHT"
    # caption: Each rectangle is as large as its share of the rows
    When user hovers over the "leaf F | RA" area of tree map viewer
    # caption: Click a rectangle to select its rows
    And user clicks on the "leaf F | RA" area of tree map viewer
    Then 1941 rows should be selected
    When the video says "A sunburst shows the same levels as rings" at the top
    And user clicks on "Add viewer" icon in toolbar
    And user types "sunburst" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be visible
    # caption: Its levels are set in the settings
    When user hovers over sunburst viewer
    And user clicks on settings icon of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    And user toggles the "DIS_POP" column in the column list of "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then the "hierarchy columns" reading of sunburst viewer should be "SEX, DIS_POP"
    # caption: Click a segment to select its rows
    When user clicks on the "segment M | Psoriasis" area of sunburst viewer
    Then 820 rows should be selected
    And no errors should have been logged
