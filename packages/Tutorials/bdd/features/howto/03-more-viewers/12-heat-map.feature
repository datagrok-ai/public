@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The heat map from the toolbox on cars: every value a
  colored cell, the map sorted by price with a double-click on its header, then the settings icon: the 0/1 columns
  hidden in the Order or Hide Columns dialog (Columns). cars is opened off camera; nothing is made on the server.

  Scenario: See a whole table at once with a heat map
    A heat map draws every value of the table as a colored cell, so rows that look alike stand out.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "A heat map colors every value: high red, low blue"
    And user clicks on heat map icon on toolbox
    Then heat map viewer should be added to the open tableview
    # caption: Hide the toolbox to give the map room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Double-click a header to sort the rows by it
    When user double-clicks on the "header price" area of heat map viewer
    Then the "sort column" reading of heat map viewer should be "price"
    When the video says "Sorted by price: the big, fast cars line up at the top" at the top
    And user hovers over the "header horsepower" area of heat map viewer
    # caption: The rest is in the settings
    And user hovers over heat map viewer
    And user clicks on settings icon of heat map viewer
    # caption: Columns: hide the ones you don't need
    And user expands "Columns" category in context panel
    And user clicks on first button in "Columns" property in context panel
    And user toggles the "diesel" column in the column list of Order or Hide Columns dialog
    And user toggles the "turbo" column in the column list of Order or Hide Columns dialog
    And user toggles the "two.doors" column in the column list of Order or Hide Columns dialog
    And user toggles the "hatchback" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then the "columns shown" reading of heat map viewer should be 13
    And no errors should have been logged
