@journey @users @realizes:views.users
Feature: The Users view
  Browse > Platform > Users as an administrator sees it: the list and its toolbar, the view modes,
  search and filters, the New dialogs up to the point where they would create someone, and what a
  user's context menu, context panel and profile show. Translated from files/TestTrack/User
  groups/users_manual_tests.md (Users-01 to 04, 06, 08, 09, 11, 13 to 17, 22) and
  playwright-public/user groups/users.test.ts.

  The user every claim is about is made by the feature, named by the time it ran: a user can never
  be deleted, so a fixed name would be there already on the second run. The users it makes stay on
  the stand.

  The search is fuzzy — a new login also brings up every login sharing its letters — and a new user
  is first in the list before any search, so a search is claimed by the counter dropping first and
  only then by the user's link: the counter keeps its old number until the result lands. An item
  outside a search is never a claim, since a long gallery renders only what has scrolled in.

  Projects, Activity, Chats and Privileges count their items and hide while the count is 0, which it
  is for a new user, so they are claimed present. The Personal pane shows the picture and how long ago the
  user joined ("just now", "a minute ago": grok_user_meta.dart), not the name, email and login the
  manual case lists — nothing there reads as a fact about this user, so the pane is claimed shown.

  Not translated: Users-05 and 07 (creating users) are users-create.feature; Users-18 to 21 are
  users-manage.feature. Users-10 (#tag search): no user carries a tag to find, and the manual case's
  own check is only that nothing fails. Users-12 (sorting): the sort menu marks its field, but
  nothing on the page says in which order the gallery then is.

  Background:
    Given user is logged in
    And a new user "opavlenko{time}v" with email "opavlenko+{time}v@datagrok.ai" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Platform > Users" tree node inside browse tree

  Scenario: The view opens from the Browse tree (Users-01)
    Then the "Users" view should be current
    And the page address should contain "/users"
    And gallery should be visible
    And gallery counter should be visible
    When user remembers the gallery counter
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The toolbar carries its controls (Users-02)
    Then the following elements should be visible:
      | New button                                          |
      | gallery search                                      |
      | "Switch to brief view" icon inside gallery toolbar  |
      | "Switch to card view" icon inside gallery toolbar   |
      | "Switch to grid view" icon inside gallery toolbar   |
      | "Sort list" icon inside gallery toolbar             |
      | "Toggle filters" icon inside gallery toolbar        |
      | "Refresh" icon inside gallery toolbar               |
    And New button should be enabled
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The view-mode icons switch the gallery's render mode (Users-11)
    Then "Switch to brief view" icon inside gallery toolbar should be selected
    And the gallery should be in brief mode
    When user clicks on "Switch to card view" icon inside gallery toolbar
    Then the gallery should be in card mode
    And "Switch to card view" icon inside gallery toolbar should be selected
    And "Switch to brief view" icon inside gallery toolbar should not be selected
    When user clicks on "Switch to grid view" icon inside gallery toolbar
    Then the gallery should be in grid mode
    And "Switch to grid view" icon inside gallery toolbar should be selected
    And "Switch to card view" icon inside gallery toolbar should not be selected
    When user clicks on "Switch to brief view" icon inside gallery toolbar
    Then the gallery should be in brief mode
    And "Switch to brief view" icon inside gallery toolbar should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Searching by login narrows the list, and clearing brings the rest back (Users-09)
    When user types "opavlenko{time}v" into gallery search
    Then the gallery counter should be lower than remembered
    And "opavlenko{time}v" link in gallery should be visible
    And the page address should contain "?q=opavlenko{time}v"
    When user remembers the gallery counter
    And user clears gallery search
    Then the gallery counter should be higher than remembered
    And the page address should not contain "?q="
    When user remembers the gallery counter
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: New offers a user, a service user and an invitation (Users-03)
    When user clicks on New button
    Then the open menu should list "User..."
    And the open menu should list "Service User..."
    And the open menu should list "Invite a Friend..."
    When user presses Escape
    Then context menu should be hidden
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The new user dialog asks for four fields, and Cancel closes it (Users-04)
    When user clicks on New button
    And user picks "User..." from the open menu
    Then "Create new user" dialog should be visible
    And the following elements should be visible:
      | Email input in "Create new user" dialog        |
      | Login input in "Create new user" dialog        |
      | "First Name" input in "Create new user" dialog |
      | "Last Name" input in "Create new user" dialog  |
    And OK button in "Create new user" dialog should be disabled
    When user clicks on CANCEL button in "Create new user" dialog
    Then the "Create new user" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The new user dialog refuses a bad email and a bad login (Users-06)
    When user clicks on New button
    And user picks "User..." from the open menu
    Then OK button in "Create new user" dialog should be disabled
    When user types "not-an-email" into Email input in "Create new user" dialog
    And user types "Opavlenko+Bad" into Login input in "Create new user" dialog
    And user types "Bad" into "First Name" input in "Create new user" dialog
    And user types "Input" into "Last Name" input in "Create new user" dialog
    Then Email input in "Create new user" dialog should be invalid
    And Login input in "Create new user" dialog should be invalid
    And OK button in "Create new user" dialog should be disabled
    When user types "opavlenko+bad{time}v@datagrok.ai" into Email input in "Create new user" dialog
    Then Email input in "Create new user" dialog should be valid
    And OK button in "Create new user" dialog should be disabled
    When user types "opavlenko-bad{time}v" into Login input in "Create new user" dialog
    Then Login input in "Create new user" dialog should be valid
    And OK button in "Create new user" dialog should be enabled
    When user clicks on CANCEL button in "Create new user" dialog
    Then the "Create new user" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Invite a Friend asks for an email (Users-08)
    When user clicks on New button
    And user picks "Invite a Friend..." from the open menu
    Then "Invite a Friend" dialog should be visible
    And Email input in "Invite a Friend" dialog should be visible
    When user clicks on CANCEL button in "Invite a Friend" dialog
    Then the "Invite a Friend" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A user's context menu (Users-14)
    When user types "opavlenko{time}v" into gallery search
    Then the gallery counter should be lower than remembered
    When user opens the context menu of "opavlenko{time}v" link in gallery
    Then the open menu should list "Details"
    And the open menu should list "Chat"
    And the open menu should list "Disable..."
    And the open menu should list "Groups..."
    And the open menu should list "Roles..."
    And the open menu should list "Copy > ID"
    And the open menu should list "Add to favorites"
    When user closes the context menu
    And user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Selecting a user fills the context panel (Users-17)
    When user types "opavlenko{time}v" into gallery search
    Then the gallery counter should be lower than remembered
    When user clicks on "opavlenko{time}v" link in gallery
    Then the context panel should show "opavlenko{time}v"
    And the following elements should be visible:
      | "Personal" accordion header in context panel    |
      | "Roles" accordion header in context panel       |
      | "Member of" accordion header in context panel   |
      | "Sticky meta" accordion header in context panel |
    And the following elements should be present:
      | "Projects" accordion header in context panel    |
      | "Activity" accordion header in context panel    |
      | "Chats" accordion header in context panel       |
      | "Privileges" accordion header in context panel  |
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Double-clicking a user opens the profile (Users-15)
    When user types "opavlenko{time}v" into gallery search
    Then the gallery counter should be lower than remembered
    When user double-clicks on "opavlenko{time}v" link in gallery
    Then the "opavlenko{time}v" view should be current
    When user closes the current view
    Then the "Users" view should be current
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Details opens the same profile (Users-16)
    When user types "opavlenko{time}v" into gallery search
    Then the gallery counter should be lower than remembered
    When user picks "Details" from the context menu of "opavlenko{time}v" link in gallery
    Then the "opavlenko{time}v" view should be current
    When user closes the current view
    Then the "Users" view should be current
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  # GROK-20905: the first time the filters open on a page, the platform logs a NullError from
  # formula_lines_mixin.dart (LineChartForFL.formulaLines) once the Joined card's histogram has drawn;
  # never on a later open. This scenario holds only that claim, read after the card is up and the panel
  # closed again; what the filters show is the next scenario's, which opens them itself, so a failure
  # there is not swallowed by the tag. Both are last: a known failure stops at its failing step.
  @known-failure
  Scenario: Opening the filters for the first time logs no error (Users-13)
    When user clicks on "Toggle filters" icon inside gallery toolbar
    Then "Joined" filter card should be visible
    When user clicks on "Toggle filters" icon inside gallery toolbar
    Then no errors should have been logged

  Scenario: The filters show a card per user property (Users-13)
    When user clicks on "Toggle filters" icon inside gallery toolbar
    Then filter panel should be visible
    And "Status" filter card should be visible
    And "Joined" filter card should be visible
    When user clicks on "Toggle filters" icon inside gallery toolbar
    Then filter panel should be hidden
    And no errors should have been logged
    And no error or warning balloon should have been shown
