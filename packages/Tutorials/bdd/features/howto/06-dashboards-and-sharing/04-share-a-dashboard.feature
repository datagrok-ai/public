@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The dashboard "Demographics overview" (demog with a scatter
  plot and a bar chart) is saved off camera through the ribbon's Save dialog, the Share offer cancelled, and found
  in Browse > Dashboards. On camera its card's Share... opens the Share dialog: the fixture user bddsecond (the
  library's sharing user, DATAGROK_SHARING_LOGIN) is typed and taken from the suggestions, the access levels are
  shown and left at View and use, Send notifications is unchecked (so no notification is sent), and OK shares it.
  A line over the card's Links... link says what it gives (the Links dialog is not opened: its address is the
  stand's, localhost), and its Sharing section lists the user (the pointer rests on the grant, so the line is read
  before the cut that follows). Then the page signs in as the sharing user, off camera: a
  stand has no second real login, so the cut is made readable by a line saying whose session it is now (the
  sidebar's User icon only says "User" on hover); their first click, Dashboards in their Browse, is on camera, and a search finds the
  dashboard (My stuff > Shared with me holds other features' leftovers on a shared stand). They open it and get
  only "Save a copy" in the Save dialog: View and use does not allow saving over the original.

  The project, with its table, views, picture and the grant, is deleted when the feature starts and ends; the
  sharing user's notifications are cleared now and at feature end, and the running account signs back in first
  thing at feature end.

  Scenario: Share a dashboard
    Give people or groups access to a dashboard, and see where they find it.

    Given user is logged in
    And simple mode is off
    And the sharing user has no notifications, now and when the feature ends
    And no project named "Demographics overview" is on the server
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE    |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And user adds a bar chart viewer with:
      | splitColumnName | RACE |
    And user clicks on Save ribbon item
    And user enters "Demographics overview" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    And the "Save project" dialog should close
    And "Share Demographics overview" dialog should be visible
    And user clicks on CANCEL button in "Share Demographics overview" dialog
    And the "Share Demographics overview" dialog should close
    And user closes all views
    And the notification balloons are closed
    And the toolbox pane is hidden
    And the context panel is open
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    And "Demographics overview" gallery card should be visible
    When the video says "A new dashboard is yours alone. Share it to let others open it"
    # caption: Right-click it, then Share...
    And user picks "Share..." from the context menu of "Demographics overview" gallery card
    Then "Share Demographics overview" dialog should be visible
    # caption: Type a user, a group or an email
    When user types "bddsecond" into "User, group, or email" input in "Share Demographics overview" dialog
    And user clicks on "bddsecond" table row in "Share Demographics overview" dialog
    # caption: Choose what they may do: view, edit, or full access
    And user clicks on share access selector
    Then privilege tree should be visible
    When user hovers over "Edit" tree node inside privilege tree
    And user presses Escape
    Then privilege tree should be hidden
    And share access selector should contain text "View and use"
    # caption: Uncheck Send notifications to share without a message
    When user unchecks "Send notifications" input in "Share Demographics overview" dialog
    And user clicks on OK button in "Share Demographics overview" dialog
    Then the "Share Demographics overview" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on "Demographics overview" gallery card
    Then the context panel should show "Demographics overview"
    # caption: Links... gives its address, for the people who have access
    When user hovers over "Links..." link in Details pane in context panel
    # caption: The Sharing section shows who can open it
    And user expands Sharing section in context panel
    Then the sharing pane should list the sharing user
    When user hovers over the text "can view and use" in context panel
    Given user closes all views
    And user signs in as the sharing user
    And the sharing user should be signed in
    And simple mode is off
    And the browse panel is open
    And the context panel is closed
    And user refreshes the browse tree
    And the notification balloons are closed
    When the video says "Now signed in as the colleague it was shared with"
    # caption: They find it in their own Browse > Dashboards
    And user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    Then "Demographics overview" gallery card should be visible
    # caption: Double-click to open it
    When user double-clicks on "Demographics overview" gallery card and the project opens
    Then the "demog" view should be current
    And scatter plot viewer should be visible
    And bar chart viewer should be visible
    # caption: With View and use they cannot save over it, only save a copy
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Name text input in "Save project" dialog should have value "Copy of Demographics overview"
    When user hovers over "Save original project" text in "Save project" dialog
    And user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged
