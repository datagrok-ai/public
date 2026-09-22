@connections @journey
Feature: A connection in the connections browser and its context panel
  A provider's "Browse connections" view lists its connections in a gallery with a search; a card
  opens the connection on the context panel — its details, sharing, activity, chat, and the menu
  behind the header's arrow. Translated from TestTrack Connections/browser.md
  (playwright-public connections/04-browser.test.ts).

  The subject is a Postgres connection saved through the API without credentials,
  BDD-Conn-Browser-{run}; the share, the chat and the connection are removed at feature end and
  checked gone (a chat goes before its connection — a chat must not outlive what it is about).

  Not translated, and why: the md's "Filter templates (magic wand)" icon — the view has no such icon
  any more (its toolbar has a plain filter). The old spec's checks that proved nothing are not
  restored: it granted the share through the API rather than the dialog, and its Activity regex
  matched the connection's own name.

  A tree row is wider than the browse panel shows, and the context-menu gesture aims at the row's
  middle — past the panel's edge for a long name, where no menu opens; the rows are right-clicked
  instead, which lands inside them (library candidate: aim the context menu at the visible part).

  Background:
    Given user is logged in
    And the browse panel is open
    And Databases tree node inside browse tree is expanded
    And a "Postgres" connection named "BDD-Conn-Browser-{run}" is on the server
    And the context panel is open

  Scenario: The connections view searches by name and opens a card on the context panel
    When user right-clicks on Databases---Postgres tree node inside browse tree
    And user picks "Browse connections" from the open menu
    Then the "Postgres" view should be current
    And the page address should contain "/connections/Postgres"
    When user types "BDD-Conn-Browser-{run}" into gallery search
    Then "BDD-Conn-Browser-{run}" link in gallery should be visible
    And there should be 1 visible link in gallery
    When user clicks on "BDD-Conn-Browser-{run}" link in gallery
    Then the context panel should show "BDD-Conn-Browser-{run}"
    And Details section in context panel should contain the text "db.datagrok.ai"
    And Details section in context panel should contain the text "northwind"
    And Details section in context panel should contain the text "54322"
    And no errors should have been logged

  Scenario: The header arrow opens the connection's menu
    When user clicks on "context-arrow-down" icon in context panel
    Then the open menu should list "Test connection"
    And the open menu should list "Share..."
    And the open menu should list "Delete..."
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Sharing the connection shows on its Sharing pane
    When user picks "Share..." from the context menu of "BDD-Conn-Browser-{run}" link in gallery
    Then "Share BDD-Conn-Browser-{run}" dialog should be visible
    And "Share BDD-Conn-Browser-{run}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDD-Conn-Browser-{run}" dialog
    And user unchecks "Send notifications" input in "Share BDD-Conn-Browser-{run}" dialog
    And user clicks on OK button in "Share BDD-Conn-Browser-{run}" dialog
    Then the "Share BDD-Conn-Browser-{run}" dialog should close
    When user clicks on "BDD-Conn-Browser-{run}" link in gallery
    Then the sharing pane should list the sharing user
    And no error or warning balloon should have been shown

  Scenario: The Activity pane records the connection
    Then Activity section in context panel should be present
    And no errors should have been logged

  Scenario: A chat message is posted on the connection and removed with it
    Given Chats section in context panel is expanded
    When user types "BDD chat {run}" into chat post input in context panel
    And user presses Enter in chat post input in context panel
    Then Chats section in context panel should contain the text "BDD chat {run}"
    And the "BDD-Conn-Browser-{run}" connection should have a chat on the server
    When user deletes the chat of the "BDD-Conn-Browser-{run}" connection
    Then the "BDD-Conn-Browser-{run}" connection should have no chat on the server
    And no errors should have been logged
