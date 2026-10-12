@journey @serial @realizes:sharing.share-dialog @realizes:sharing.browse-shared-with-me
Feature: A folder in your files, shared
  A single file has no share command; a folder has Share folder.... The first share of a folder
  makes a connection for that folder, named "My files/<folder>", and shares it; the second account
  then finds the folder under My stuff > Shared with me. Translated from the TestTrack case
  Sharing/sharing-folder.

  The folder ShareFolder1 in Browse > Files > My files holds one file, data1.csv (the project's
  fixtures/cars-small.csv), written through the API; the file goes when the feature ends. The
  connection the share makes (MyFilesShareFolder1 on the server) is removed through its Delete...
  in the last scenario, and swept by name when the feature starts and ends, its grant with it; the
  folder is removed through its Delete... in the last scenario. The second account is the library's
  sharing user; it signs in on the feature's page, and the share is made with Send notifications off.

  Kept without: that the reopened Share folder dialog lists the second account (no step reads a
  grant row of the sharing user in a dialog). What the second account may do with the files (read,
  not write) and that files outside the folder stay closed to it is an API test.

  Known failure: GitHub #4268 "Sharing: My files: error on opening the shared folder" — the second
  account gets a "Credentials required" view instead of the file list (localhost, 2026-10-12). The
  folder's view is first opened from Shared with me and claimed current in an ordinary scenario, and
  closed; the known-failure scenario opens it again the same way and claims only its content.
  Both opens go through the node under My stuff > Shared with me: the second account also sees the
  folder's connection under Files.

  Serial: it signs in as the second account on its page.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And the user's own connection "MyFilesShareFolder1" is removed now and at feature end
    And the "fixtures/cars-small.csv" file of the project is in the user's files as "ShareFolder1/data1.csv"
    And the browse panel is open
    When user refreshes the browse tree
    Given Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded

  Scenario: Share folder opens and cancels cleanly (GROK-19264)
    When user picks "Share folder..." from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree
    Then "Share folder" dialog should be visible
    And link icon in "Share folder" dialog should be hidden
    And "User, group, or email" input in "Share folder" dialog should be visible
    And share access selector should contain text "View and use"
    And OK button in "Share folder" dialog should be visible
    And CANCEL button in "Share folder" dialog should be visible
    When user clicks on CANCEL button in "Share folder" dialog
    Then the "Share folder" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown
    And 0 connections named "MyFilesShareFolder1" should be on the server
    When user refreshes the browse tree
    Given Files---My-files tree node inside browse tree is expanded
    Then Files---My-files---ShareFolder1 tree node inside browse tree should be visible
    And "My files/ShareFolder1" tree node inside browse tree should be absent

  Scenario: The folder is shared with the second account (GROK-12795, GROK-11266)
    When user picks "Share folder..." from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree
    Then "Share folder" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share folder" dialog
    And user unchecks "Send notifications" input in "Share folder" dialog
    And user clicks on OK button in "Share folder" dialog
    Then the "Share folder" dialog should close
    And an info balloon containing "Shared" should have been shown
    And no error or warning balloon should have been shown
    When user refreshes the browse tree
    Then "My files/ShareFolder1" tree node inside browse tree should be visible
    When user clicks on "My files/ShareFolder1" tree node inside browse tree
    Then data1.csv link in gallery should be visible
    Given Files---My-files tree node inside browse tree is expanded
    When user picks "Share folder..." from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree
    Then "Share folder" dialog should be visible
    And link icon in "Share folder" dialog should be visible
    When user clicks on CANCEL button in "Share folder" dialog
    Then the "Share folder" dialog should close

  Scenario: The second account finds the folder under Shared with me
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    And user double-clicks on "My files/ShareFolder1" shared by the running account
    Then the "My files/ShareFolder1" view should be current
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current

  @known-failure
  Scenario: The second account sees the folder's files, not a credentials form (#4268)
    When user double-clicks on "My files/ShareFolder1" shared by the running account
    Then data1.csv link in gallery should be visible
    And "Credentials required" text should be absent

  Scenario: The owner deletes the folder's connection and the folder
    When user signs in as themselves again
    Then the running account should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given Files tree node inside browse tree is expanded
    When user picks "Delete..." from the context menu of "My files/ShareFolder1" tree node inside browse tree
    Then "Are you sure?" dialog should contain text 'Delete connection "My files/ShareFolder1"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 connections named "MyFilesShareFolder1" should be on the server
    Given Files---My-files tree node inside browse tree is expanded
    When user picks "Delete..." from the context menu of Files---My-files---ShareFolder1 tree node inside browse tree
    Then "Are you sure?" dialog should contain text "ShareFolder1"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user refreshes the browse tree
    Given Files---My-files tree node inside browse tree is expanded
    Then "My files/ShareFolder1" tree node inside browse tree should be absent
    And Files---My-files---ShareFolder1 tree node inside browse tree should be absent
