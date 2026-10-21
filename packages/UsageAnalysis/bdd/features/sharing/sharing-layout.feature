@journey @serial @realizes:sharing.share-dialog @realizes:views.layouts
Feature: A layout, shared and applied by the other user
  A layout the owner shares from its card in the Layouts pane, seen from the other side: the second
  account opening a table of the same columns gets the layout's card and applies it with a click.
  Translated from playwright-public/Sharing share-layout-permissions (block D, whose second-account
  steps were API reads); the owner's Share dialog on a layout card is layouts/layout-names-tags, and
  the dialog's own controls and the Sharing pane before a share (blocks A and B) are sharing-dialog
  (on a query); the Advanced editor of block C is gone from the dialog (GROK-20322).

  The layout is this feature's own, ShareLayout-{time}, saved with Save in the Layouts pane from
  demog's first 100 rows with a scatter plot of WEIGHT by HEIGHT (a layout saved through the JS API
  is not offered in the pane), renamed in the context panel, and deleted when the feature ends. The
  second account is the library's sharing user; it signs in on the feature's page, and the share is
  made with Send notifications off.

  Not translated, and why: that the second account may not save over or delete the layout, and loses
  it after the revoke, is the server's answer (ApiTests, Dapi: permissions: "a layout", "a revoke").

  Serial: it signs in as the second account on its page.

  Background:
    Given user is logged in
    And the sharing user has no notifications, now and when the feature ends
    And the layouts named "ShareLayout-{time}" are deleted when the feature ends
    And the layouts named "share-layout-{time}" are deleted when the feature ends
    And user opens demog dataset keeping the first 100 rows as "share-layout-{time}"
    And user adds a scatter plot viewer with:
      | xColumnName | HEIGHT |
      | yColumnName | WEIGHT |
    And the toolbox pane is shown
    And the context panel is open
    And Layouts accordion header in toolbox is expanded

  Scenario: The owner saves a layout from the Layouts pane and shares it from its card
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    When user clicks on Save button in layouts pane
    Then the context panel should show "share-layout-{time}"
    When user types "ShareLayout-{time}" into Name details field in context panel
    And user presses Enter
    Then "ShareLayout-{time}" layout card should be visible
    When user picks "Share..." from the context menu of "ShareLayout-{time}" layout card
    And user picks the sharing user in "User, group, or email" input in "Share ShareLayout-{time}" dialog
    And user unchecks "Send notifications" input in "Share ShareLayout-{time}" dialog
    And user clicks on OK button in "Share ShareLayout-{time}" dialog
    Then the "Share ShareLayout-{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account gets the layout's card for a table of the same columns and applies it
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given user opens demog dataset keeping the first 100 rows as "share-layout-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    Then "ShareLayout-{time}" layout card should be visible
    And the open tableview should have 0 scatter plot viewers
    When user clicks on "ShareLayout-{time}" layout card
    Then the open tableview should have 1 scatter plot viewer
    And "xColumnName" property of scatter plot viewer should be "HEIGHT"
    And "yColumnName" property of scatter plot viewer should be "WEIGHT"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user signs in as themselves again
    Then the running account should be signed in
