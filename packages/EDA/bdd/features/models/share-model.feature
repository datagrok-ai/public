@journey @serial @eda @realizes:sharing.share-dialog @realizes:sharing.browse-shared-with-me
Feature: Sharing a predictive model
  A model its owner trains and saves, shared with a second account and taken back. Translated from
  playwright-public/Sharing/share-model-permissions.md and the package's
  playwright/share-model-permissions.test.ts: the owner's side through the Share dialog, claimed in the
  Sharing pane of the context panel, and the second account's side (blocks D.2 and D.3) — it finds
  the model under My stuff > Shared with me and applies it to its own demog-1000. The old spec's
  second-account steps were API reads.

  The model is the case's own example, SEX predicted by HEIGHT and WEIGHT, trained on demog-1000. HEIGHT
  has missing values, so the view offers an engine only once Ignore missing is on. The picker lists the
  table in its column order, HEIGHT sixth and WEIGHT seventh.

  Not translated, and why: that the second account is refused edit, delete and re-share (block F) and
  loses the model after the revoke (block G) is the server's answer (ApiTests, Dapi: permissions,
  "a model"), and so is block E (no dependent entity shared along); the Advanced editor link of block
  C is no longer in the dialog (GROK-20322). The second account is DATAGROK_SHARING_LOGIN or the
  "bddsecond" user the library's setup creates, and must not be an administrator, who sees the model
  unshared; it signs in on the feature's page.

  Serial: it signs in as the second account on its page.

  The Share dialog of a model is ready only once it has fetched the model's project, and says so by
  listing the owner's grant ("Full access"); an OK pressed before that fails with "Not initialized"
  and leaves the dialog open. A space is a project itself, so its dialog is ready at once. The context
  panel opens only after the model is saved: while the Train Model view is current the panel holds
  the view's own properties.

  Background:
    Given user is logged in
    And no predictive model named "BDD-Share-Model-{run}" is on the server

  Scenario: The owner trains a model and saves it
    Given user opens demog-1000 dataset
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "SEX" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    Then the "text of cell 6 of __name" reading of grid viewer in "Select columns..." dialog should be "HEIGHT"
    And the "text of cell 7 of __name" reading of grid viewer in "Select columns..." dialog should be "WEIGHT"
    When user clicks on the "cell 6 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 7 of x" area of grid viewer in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "2 checked"
    When user clicks on OK button in "Select columns..." dialog
    And user checks "Ignore missing" input
    And user selects "Eda: XGBoost" in "Model Engine" input
    Then model preview should be ready
    And "Accuracy" table row should be visible
    When user clicks on Save button
    And user enters "BDD-Share-Model-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then 1 predictive model named "BDD-Share-Model-{run}" should be on the server

  Scenario: The saved model is the owner's alone
    Given the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on "Predictive models" tree node inside browse tree
    Then the "Models" view should be current
    When user clicks on "BDD-Share-Model-{run}" label in gallery
    Then the context panel should show "BDD-Share-Model-{run}"
    And the sharing pane should not list the sharing user

  Scenario: The Share dialog asks who, how much, and whether to notify
    When user picks "Share..." from the context menu of "BDD-Share-Model-{run}" label in gallery
    Then "Share BDD-Share-Model-{run}" dialog should be visible
    And "User, group, or email" input in "Share BDD-Share-Model-{run}" dialog should be visible
    And share access selector should contain text "View and use"
    And "Send notifications" input in "Share BDD-Share-Model-{run}" dialog should be hidden
    And "Share BDD-Share-Model-{run}" dialog should not contain text "will also be shared"
    When user clicks on CANCEL button in "Share BDD-Share-Model-{run}" dialog
    Then "Share BDD-Share-Model-{run}" dialog should be hidden
    When user clicks on "BDD-Share-Model-{run}" label in gallery
    Then the sharing pane should not list the sharing user

  Scenario: The model is shared with the second account to view and use
    When user picks "Share..." from the context menu of "BDD-Share-Model-{run}" label in gallery
    Then "Share BDD-Share-Model-{run}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDD-Share-Model-{run}" dialog
    Then "Send notifications" input in "Share BDD-Share-Model-{run}" dialog should be visible
    And "Share BDD-Share-Model-{run}" dialog should not contain text "will also be shared"
    When user unchecks "Send notifications" input in "Share BDD-Share-Model-{run}" dialog
    And user clicks on OK button in "Share BDD-Share-Model-{run}" dialog
    Then the "Share BDD-Share-Model-{run}" dialog should close
    And no error or warning balloon should have been shown
    When user clicks on "BDD-Share-Model-{run}" label in gallery
    Then the sharing pane should list the sharing user

  Scenario: The second account finds the model under Shared with me and applies it
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user refreshes the browse tree
    Given "My stuff" tree node inside browse tree is expanded
    And "My stuff > Shared with me" tree node inside browse tree is expanded
    When user expands "." shared by the running account
    Then "BDD-Share-Model-{run}" shared by the running account should be visible
    Given user opens demog-1000 dataset
    When user picks "ML > Models > Apply Model..." from the top menu
    Then "Apply predictive model" dialog should be visible
    When user selects the predictive model "BDD-Share-Model-{run}" in Model input in "Apply predictive model" dialog
    And user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And 1 new column should have been added
    And a new column matching "^SEX" should have been added
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user signs in as themselves again
    Then the running account should be signed in
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on "Predictive models" tree node inside browse tree
    Then the "Models" view should be current

  Scenario: The owner takes the share back
    When user picks "Share..." from the context menu of "BDD-Share-Model-{run}" label in gallery
    Then the access level of the sharing user in "Share BDD-Share-Model-{run}" dialog should be "View and use"
    When user removes the sharing user from "Share BDD-Share-Model-{run}" dialog
    And user clicks on OK button in "Share BDD-Share-Model-{run}" dialog
    Then the "Share BDD-Share-Model-{run}" dialog should close
    And no error or warning balloon should have been shown
    When user clicks on "BDD-Share-Model-{run}" label in gallery
    Then the sharing pane should not list the sharing user
    And no errors should have been logged
