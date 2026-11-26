@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The space "Clinical team" is made through the API, and the
  dashboard "Demographics overview" (demog with a scatter plot and a bar chart) is saved off camera through the
  ribbon's Save dialog, the Share offer cancelled, and found in Browse > Dashboards. On camera the card's Move to
  Space... moves the dashboard into the space, Browse > Spaces is expanded and the space opened with a double click,
  its view lists the dashboard, the space is shared from its node in
  Browse > Spaces with the fixture user bddsecond (the library's sharing user) at View and use with Send
  notifications off, and the dashboard's own Share dialog shows the access as inherited from the space.

  The project and the space, the space's grant with it, are deleted when the feature starts and ends; the sharing
  user's notifications are cleared now and at feature end.

  Scenario: Share a space with your team
    Put dashboards in a space and share the space once: everything in it goes to the same people.

    Given user is logged in
    And simple mode is off
    And the sharing user has no notifications, now and when the feature ends
    And no project named "Demographics overview" is on the server
    And no space named "Clinical team" is on the server
    And a space named "Clinical team" is on the server
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
    And the context panel is closed
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    And "Demographics overview" gallery card should be visible
    When the video says "To share many dashboards at once, put them in a space"
    # caption: Move the dashboard into your team's space
    And user picks "Move to Space..." from the context menu of "Demographics overview" gallery card
    Then "Move to space" dialog should be visible
    When user selects "Clinical team" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    # caption: A space keeps a team's dashboards, files and queries together
    When user expands Spaces tree node inside browse tree
    And user double-clicks on "Spaces > Clinical team" tree node inside browse tree
    Then the "Clinical team" view should be current
    And "Demographics overview" link in gallery should be visible
    # caption: Share the space once
    When user picks "Share..." from the context menu of "Spaces > Clinical team" tree node inside browse tree
    Then "Share Clinical team" dialog should be visible
    When user types "bddsecond" into "User, group, or email" input in "Share Clinical team" dialog
    And user clicks on "bddsecond" table row in "Share Clinical team" dialog
    And user unchecks "Send notifications" input in "Share Clinical team" dialog
    And user clicks on OK button in "Share Clinical team" dialog
    Then the "Share Clinical team" dialog should close
    # caption: Now look at the dashboard inside
    When user picks "Share..." from the context menu of "Demographics overview" link in gallery
    Then "Share Demographics overview" dialog should be visible
    And "Share Demographics overview" dialog should contain text "Inherited from"
    And "Share Demographics overview" dialog should contain text "Clinical team"
    # caption: It has the same access, inherited from the space
    When user hovers over "Share Demographics overview" dialog
    Then no errors should have been logged
