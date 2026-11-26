@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The Tile Viewer from the Add viewer gallery on demog
  (dragged by its title bar to the bottom on camera, so the lanes have room): a card per row, then the settings icon: the cards split
  into lanes by SEVERITY, and the first row's card dragged from None to High, which writes High into its cell
  (a change to the open table only). Then the context panel hidden with the status bar's
  toggle and the Form viewer from the gallery: a click on a card (row 2's, first in the None lane; row 3's sat
  half hidden under row 1's in the High lane, which scrolled under the press and left the outline on row 1's card)
  shows that row in the form, and its next-row arrow steps to the next one. demog is opened off camera; nothing is made on the
  server.

  Scenario: Show rows as cards with the Tile Viewer and the Form
    The Tile Viewer shows every row as a card and the Form shows one row at a time; both edit the table.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is hidden
    When the video says "The Tile Viewer shows every row as a card"
    And user clicks on "Add viewer" icon in toolbar
    And user types "tile" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Tile Viewer" card in "Add Viewer" dialog
    Then Tile Viewer viewer should be visible
    # caption: Drag its title bar to the bottom, so the lanes have room
    When user docks Tile Viewer viewer to the bottom edge of the view
    Then Tile Viewer viewer should be docked along the bottom edge of the view
    # caption: Split the cards into lanes by a column
    When user hovers over Tile Viewer viewer
    And user clicks on settings icon of Tile Viewer viewer
    And user selects "SEVERITY" in "Lanes" property in context panel
    Then the "lanes" reading of Tile Viewer viewer should be 5
    # caption: Drag a card to another lane to change its value
    When user drags the card of row 1 of Tile Viewer viewer to the "High" lane
    Then the value of "SEVERITY" column in row 1 should be "High"
    # caption: Hide the context panel for more room
    When user clicks on "Context Panel" status bar toggle
    Then context panel should be hidden
    When the video says "The Form viewer shows one row at a time" at the top
    And user clicks on "Add viewer" icon in toolbar
    And user types "form" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Form" card in "Add Viewer" dialog
    Then form viewer should be visible
    # caption: Click a card: the form shows its row
    When user clicks on the "tile of row 2" area of Tile Viewer viewer
    Then the "row" reading of form viewer should be 2
    # caption: The arrows step through the rows
    When user clicks on the "next row" area of form viewer
    Then the "row" reading of form viewer should be 3
    And no errors should have been logged
