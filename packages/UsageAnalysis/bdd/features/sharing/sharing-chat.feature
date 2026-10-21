@journey @serial @realizes:sharing.share-dialog @realizes:views.queries
Feature: Two people in the chat of a shared entity
  The owner posts in the Chats pane of an entity and shares it; the second account sees the post in
  the same pane and replies, and the owner sees the reply. Translated from the TestTrack case
  Chat/projects-chat-collaboration (steps 1-4, 6 and 7) and playwright-public/Chat
  projects-chat-collaboration.

  The md chats on a project; this feature chats on a query, ShareChat1 on System:Datagrok, since the
  library deletes a query with its chats when the feature starts and ends, and has no step that removes
  a project's chat (a chat must not outlive what it is about). The Chats pane is the same for both.
  The second account is the library's sharing user; it signs in on the feature's page, sees the query
  under My stuff > Shared with me (which lists only what was shared with it), opens it from the
  Queries view as the owner does, and the share is made with Send notifications off.

  Not translated, and why: the sharing notification (step 3, 5) is claimed by PowerPack
  home/home-as-sharing-user; an empty chat before the first post (step 1), who wrote a message, its
  avatar and time, the order of the messages (steps 2, 6, 7), mentions
  (8), editing and deleting a message (9, 10) need readings no step has (sharing/MISSING.md); that the
  history is stored and read back across sessions (11) is the datlas test dapi/chats_test.dart.

  Serial: it signs in as the second account on its page, and the Queries view's search text is the
  account's own setting.

  Background:
    Given user is logged in
    And a query "ShareChat1" on "System:Datagrok" reads "select 1 as x"
    And the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareChat1" into gallery search
    Then "ShareChat1" gallery card should be visible

  Scenario: The owner posts in the query's chat and shares the query
    When user clicks on "ShareChat1" gallery card
    Then the context panel should show "ShareChat1"
    Given Chats section in context panel is expanded
    When user types "First comment from A {run}" into chat post input in context panel
    And user presses Enter in chat post input in context panel
    Then Chats section in context panel should contain the text "First comment from A {run}"
    When user picks "Share..." from the context menu of "ShareChat1" gallery card
    Then "Share ShareChat1" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share ShareChat1" dialog
    And user unchecks "Send notifications" input in "Share ShareChat1" dialog
    And user clicks on OK button in "Share ShareChat1" dialog
    Then the "Share ShareChat1" dialog should close
    And an info balloon containing "Shared" should have been shown
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account sees the post and replies
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the context panel is open
    And the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "ShareChat1" shared by the running account should be visible
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareChat1" into gallery search
    And user clicks on "ShareChat1" gallery card
    Then the context panel should show "ShareChat1"
    Given Chats section in context panel is expanded
    Then Chats section in context panel should contain the text "First comment from A {run}"
    When user types "Reply from B {run}" into chat post input in context panel
    And user presses Enter in chat post input in context panel
    Then Chats section in context panel should contain the text "Reply from B {run}"
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The owner sees the reply
    When user signs in as themselves again
    Then the running account should be signed in
    Given the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Queries tree node inside browse tree
    Then the "Queries" view should be current
    When user enters "ShareChat1" into gallery search
    And user clicks on "ShareChat1" gallery card
    Then the context panel should show "ShareChat1"
    Given Chats section in context panel is expanded
    Then Chats section in context panel should contain the text "First comment from A {run}"
    And Chats section in context panel should contain the text "Reply from B {run}"
    When user clears gallery search
    Then gallery search should have value ""
    And no errors should have been logged
