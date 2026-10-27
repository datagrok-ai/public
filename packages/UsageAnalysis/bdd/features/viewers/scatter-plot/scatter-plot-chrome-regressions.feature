@viewers @realizes:viewers.chrome
Feature: Viewer chrome regressions — the context panel, the title and layouts
  Closed issues about what every viewer inherits, claimed on the scatter plot (the viewer the issues
  were reported on): the gear of a viewer on the first table shows its own properties after a viewer
  on the second table was edited (public#667); after a property is edited in the context panel, the
  gear of another viewer switches the panel to that viewer (public#719); a viewer stacked over the grid
  and then closed does not come back as a floating window from the layout (public#1745); a title typed
  into the viewer's header becomes its Title and comes back from a layout (public#2925); a layout with a
  viewer bound to another table applies over the tables opened anew with no warning, and the viewer
  stays bound to that table (public#338). A stacked viewer has no title bar, so it is closed with
  General > Close from its menu. A title bar shows no text while the title is empty, so the title is
  set to "Plot" before it is double-clicked. Each scenario starts from demog-1000 and adds its own
  viewers; the layouts saved to the server go at feature end.
  Save as PNG (public#942) and Control+C (public#3052) make the picture but log three errors for a
  font the page cannot fetch: see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And the context panel is open

  # public#667
  @realizes:public#667
  Scenario: The gear of a viewer on the first table shows its properties after a viewer on the second table was edited (public#667)
    Given user adds a scatter plot viewer
    And user opens cars dataset
    And user adds a histogram viewer
    When user clicks on settings icon of histogram viewer
    Then the context panel should show "Histogram"
    Given "X Axis" category in context panel is expanded
    When user enters "30" in "Bins" property in context panel
    Then "Bins" property of histogram viewer should be "30"
    Given user switches to the "demog-1000" table view
    When user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    And "Bins" property in context panel should be absent
    And no errors should have been logged

  # public#719
  @realizes:public#719
  Scenario: After a property is edited in the context panel, the gear of another viewer switches the panel to it (public#719)
    Given user adds a scatter plot viewer
    And user adds a histogram viewer
    When user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    Given "Color" category in context panel is expanded
    When user selects "RACE" in "Color" property in context panel
    Then "Color" property of scatter plot viewer should be "RACE"
    When user clicks on settings icon of histogram viewer
    Then the context panel should show "Histogram"
    And "Marker Type" property in context panel should be absent
    And no errors should have been logged

  # public#1745
  @realizes:public#1745
  Scenario: A viewer stacked over the grid and closed does not come back as a floating window from the layout (public#1745)
    Given user adds a scatter plot viewer
    And user adds a histogram viewer
    When user stacks scatter plot viewer onto grid as a tab
    Then the viewers "demog-1000, Scatter plot" should be the tabs of one panel
    When user picks "General > Close" from the viewer menu of scatter plot viewer
    Then the open tableview should have 0 scatter plot viewers
    And the open tableview should have 1 histogram viewer
    When user saves the layout of the current table view to the server
    And user picks "General > Close" from the viewer menu of histogram viewer
    Then the open tableview should have 0 histogram viewers
    When user loads the saved layout
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 0 scatter plot viewers
    And grid should be visible
    And no errors should have been logged

  # public#2925
  @realizes:public#2925
  Scenario: A title typed into the viewer's header becomes its Title and comes back from a layout (public#2925)
    Given user adds a scatter plot viewer
    And user sets "Title" property of scatter plot viewer to "Plot"
    Then title of scatter plot viewer should have text "Plot"
    When user double-clicks on title of scatter plot viewer
    And user types "Heights" at the caret
    And user presses Enter
    Then title of scatter plot viewer should have text "Heights"
    And "Title" property of scatter plot viewer should be "Heights"
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user loads the saved layout
    Then "Title" property of scatter plot viewer should be "Heights"
    And title of scatter plot viewer should have text "Heights"
    And no errors should have been logged

  # public#338
  @realizes:public#338
  Scenario: A layout with a viewer bound to another table applies over the tables opened anew with no warning (public#338)
    Given user opens spgi-100 dataset
    And user adds a scatter plot viewer
    When user clicks on settings icon of scatter plot viewer
    And user selects "demog-1000" in "Table" property in context panel
    Then scatter plot viewer should be bound to table "demog-1000"
    When user saves the layout of the current table view to the server
    And user closes all views
    And user opens demog-1000 dataset
    And user opens spgi-100 dataset
    And user loads the saved layout
    Then scatter plot viewer should be bound to table "demog-1000"
    And no error or warning balloon should have been shown
    And no errors should have been logged
