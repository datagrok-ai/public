@guide @tutorial
Feature: How-to · Dashboards
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). demog with a scatter plot, a histogram and a bar chart (arranged off
  camera) is saved through the ribbon's SAVE as the dashboard "Demographics overview"; the Share dialog the first save
  opens is cancelled, which keeps it private. Every move after that is filmed: the view's tab is closed (which closes
  the project, and the platform shows Browse again by itself: a click on the sidebar's Browse icon there would close
  it), Browse > Dashboards opens the gallery, a search finds the card, the status bar's Context Panel toggle shows its
  details, My stuff > Recent lists it too. Then (it absorbed the former "Keep a dashboard in Favorites" clip) the
  card's Add To Favorites > Only for me stars it (the menu stays open after the pick, a checkable item: Escape closes
  it, on camera), and the sidebar's Favorites panel lists it and opens it with a double click.

  The project, with its table, views and picture, is deleted when the feature starts and ends; the favorite is taken
  out now and at feature end ("is not in favorites" puts the account's favorites back as they were).

  Scenario: Save a dashboard and find it again
    Save a table with its viewers as a dashboard, then find it again in Browse and in your Favorites.

    Given user is logged in
    And simple mode is off
    And no project named "Demographics overview" is on the server
    And "Demographics overview" is not in favorites
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE    |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And user adds a histogram viewer with:
      | valueColumnName | HEIGHT |
    And user adds a bar chart viewer with:
      | splitColumnName | RACE |
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "A dashboard is a table and its viewers, saved together"
    # caption: SAVE puts the table and its viewers on the server
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # caption: Give it a name you will recognize
    When user enters "Demographics overview" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "Demographics overview" should be on the server
    And "Share Demographics overview" dialog should be visible
    # caption: Next it offers to share. Until you do, only you can see it
    When user clicks on CANCEL button in "Share Demographics overview" dialog
    Then the "Share Demographics overview" dialog should close
    # caption: Close it for now. It is kept on the server
    When user closes "demog" view tab
    Then no table should be left in the workspace
    And browse panel should be visible
    When the video says "Later, find it again in Browse"
    # caption: Dashboards lists every dashboard you can open. Search by name
    When user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    Then the "Projects" view should be current
    And "Demographics overview" gallery card should be visible
    # caption: The Context Panel shows details of what you click
    When user clicks on "Context Panel" status bar toggle
    Then context panel should be visible
    # caption: Click the card: who made it and when
    When user clicks on "Demographics overview" gallery card
    Then the context panel should show "Demographics overview"
    # caption: What you used lately is also under My stuff > Recent
    When user expands "My stuff" tree node inside browse tree
    And user expands "My stuff > Recent" tree node inside browse tree
    Then "My stuff > Recent > Demographics overview" tree node inside browse tree should be visible
    # caption: Keep the ones you open often in Favorites: only for you, or for a group
    When user picks "Add To Favorites > Only for me" from the context menu of "Demographics overview" gallery card
    Then "Demographics overview" should be in favorites on the server
    # caption: The menu stays open for more groups. Escape closes it
    When user presses Escape
    # caption: The star on the left opens your favorites
    When user clicks on favorites sidebar tab
    Then favorites panel should be visible
    And "Demographics overview" tree node in favorites panel should be visible
    # caption: Double-click to open it, with all its viewers
    When user double-clicks on "Demographics overview" tree node in favorites panel and the project opens
    Then the "demog" view should be current
    And scatter plot viewer should be visible
    And histogram viewer should be visible
    And bar chart viewer should be visible
    And no errors should have been logged
