@journey
Feature: The table view's keyboard shortcuts open, clear and save what the help says
  The shortcuts of `help/datagrok/navigation/shortcuts.md` (and the Console and Variables sections of
  `navigation.md`) that no other feature presses, each claimed by what the user sees open or change: Alt+X
  hides and shows the Toolbox, Alt+V docks and closes the Variables panel, the tilde opens the Console, Alt+S
  opens Settings, Ctrl+Shift+V opens a second view of the same table, Ctrl+Shift+R clears the view down to
  the grid (and Ctrl+Z gives the viewers back), Ctrl+L docks the layout and viewer suggestion panels, and
  Ctrl+S saves the view's layout to the server, named after the table, and shows it in the context panel.

  Left out on purpose, one line each: F4, F7, Alt+T, Ctrl+O, Ctrl+H and Ctrl+Z are pressed by their own
  features; F1 toggling the help panel is `biostructure-viewer/file-open-and-preview.feature`'s, and the help
  following the current object needs a stand that serves the help pages (MISSING.md); F11 and Alt+F ask the
  browser for full screen, which no step reads and the shared page would keep for the next features (MISSING.md);
  Alt+C is `column-manager.feature`'s; which layouts the suggestion panel offers is
  `layouts/layout-edge-cases.feature`'s, so Ctrl+L claims only the panels it docks; the "Ctrl+Z to undo" hint
  Ctrl+Shift+R shows is a top-right info label, not a balloon, and no step reads it (MISSING.md).

  One journey on a three-row table written here, named with the run's time so the layout Ctrl+S saves
  is this run's; that layout is deleted when the feature ends. Every scenario puts back what it opened,
  except the two Ctrl+L panels, which no step closes (MISSING.md, needs a step).

  Background:
    Given user is logged in
    And the layouts named "bdd-keys-{time}" are deleted when the feature ends
    And user opens a table "bdd-keys-{time}" with:
      | name | num |
      | a    | 1   |
      | b    | 2   |
      | c    | 3   |
    And the context panel is open
    And the toolbox pane is shown

  Scenario: Alt+X hides the Toolbox and shows it again
    Then toolbox should be visible
    When user presses Alt+X
    Then toolbox should be hidden
    And "Toolbox" status bar toggle should not be selected
    When user presses Alt+X
    Then toolbox should be visible
    And "Toolbox" status bar toggle should be selected
    And no errors should have been logged

  Scenario: Alt+V docks the Variables panel and closes it again
    Then "Variables" dock panel should be absent
    When user presses Alt+V
    Then "Variables" dock panel should be visible
    And "Variables" status bar toggle should be selected
    When user presses Alt+V
    Then "Variables" dock panel should be absent
    And no errors should have been logged

  Scenario: The tilde opens the Console with its command line
    Then console should be absent
    When user presses Backquote
    Then console should be visible
    And "Console" status bar toggle should be selected
    When user closes the console
    Then console should be hidden
    And no errors should have been logged

  Scenario: Alt+S opens the Settings
    When user presses Alt+S
    Then the current view should be a settings view
    And the "Settings - General" view should be current
    When user closes the current view
    Then the "bdd-keys-{time}" view should be current
    And no errors should have been logged

  Scenario: Ctrl+Shift+V opens a second view on the same table
    When user presses Control+Shift+V
    Then the "bdd-keys-{time} copy" view should be current
    And the open table views should be exactly "bdd-keys-{time}, bdd-keys-{time}"
    When user selects rows where "name" is "b"
    And user switches to the "bdd-keys-{time}" view
    Then only rows where "name" is "b" should be selected
    When user switches to the "bdd-keys-{time} copy" view
    And user closes the current view
    Then the open table views should be exactly "bdd-keys-{time}"
    And no errors should have been logged

  Scenario: Ctrl+Shift+R leaves only the grid, and Ctrl+Z gives the viewers back
    Given user adds a histogram viewer
    And user adds a scatter plot viewer
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 1 scatter plot viewer
    When user clicks on the "cell 1 of name" area of grid
    And user presses Control+Shift+R
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    And grid should be visible
    When user presses Control+Z
    Then the open tableview should have 1 histogram viewer
    And the open tableview should have 1 scatter plot viewer
    When user presses Control+Shift+R
    Then the open tableview should have 0 histogram viewers
    And the open tableview should have 0 scatter plot viewers
    And no errors should have been logged

  Scenario: Ctrl+S saves the layout named after the table and shows it in the context panel
    When user clicks on the "header num" area of grid
    Then the context panel should show "num"
    When user presses Control+S
    Then the context panel should show "bdd-keys-{time}"
    And context panel should contain text "Layout saved"
    Given Layouts accordion header in toolbox is expanded
    And "bdd-keys-{time}" layout card should be visible
    And no errors should have been logged

  Scenario: Ctrl+L docks the layout and viewer suggestion panels
    Then "Layout suggestions" dock panel should be absent
    When user presses Control+L
    Then "Layout suggestions" dock panel should be visible
    And "Viewer suggestions" dock panel should be visible
    And no errors should have been logged
