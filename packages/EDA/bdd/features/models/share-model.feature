@journey @eda @realizes:sharing.share-dialog
Feature: Sharing a predictive model
  A model its owner trains and saves, shared with a second account and taken back. Translated from
  playwright-public/Sharing/share-model-permissions.md and the package's
  playwright/share-model-permissions.test.ts, the way the Spaces features share a space: the owner's
  side, through the Share dialog, claimed in the Sharing pane of the context panel.

  The model is the case's own example, SEX predicted by HEIGHT and WEIGHT on demog. Both columns have
  missing values, so the view offers an engine only once Ignore missing is on. The picker lists demog
  in its column order, HEIGHT sixth and WEIGHT seventh.

  Not translated: the recipient's side of the case (blocks D.2, D.3, E, F and the recipient's half of
  G — seeing the model under Shared with me, applying it, being refused edit, delete and re-share,
  losing it after the revoke). They need a second signed-in session, and a feature has one page. The
  Advanced editor link of block C is not in the dialog on dev (2026-09-11). The second account is
  DATAGROK_SHARING_LOGIN or the "bddsecond" user the library's setup creates, as for Spaces.

  The Share dialog of a model is ready only once it has fetched the model's project, and says so by
  listing the owner's grant ("Full access"); an OK pressed before that fails with "Not initialized"
  and leaves the dialog open. A space is a project itself, so its dialog is ready at once. The context
  panel opens only after the model is saved: while the Train Model view is current the panel holds
  the view's own properties.

  Background:
    Given user is logged in
    And no predictive model named "BDD-Share-Model" is on the server

  Scenario: The owner trains a model and saves it
    Given user opens demog dataset
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
    Then "Accuracy" table row should be visible
    When user clicks on Save button
    And user enters "BDD-Share-Model" into Name input in dialog
    And user clicks on OK button in dialog
    Then 1 predictive model named "BDD-Share-Model" should be on the server

  Scenario: The saved model is the owner's alone
    Given the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on "Predictive models" tree node inside browse tree
    Then the "Models" view should be current
    When user clicks on "BDD-Share-Model" label in gallery
    Then the context panel should show "BDD-Share-Model"
    And the sharing pane should not list the sharing user

  Scenario: The Share dialog asks who, how much, and whether to notify
    When user picks "Share..." from the context menu of "BDD-Share-Model" label in gallery
    Then "Share BDD-Share-Model" dialog should be visible
    And "Share BDD-Share-Model" dialog should contain text "Full access"
    And "User, group, or email" input in "Share BDD-Share-Model" dialog should be visible
    And share access selector should contain text "View and use"
    And "Send notifications" input in "Share BDD-Share-Model" dialog should be hidden
    And "Share BDD-Share-Model" dialog should not contain text "will also be shared"
    When user clicks on CANCEL button in "Share BDD-Share-Model" dialog
    Then "Share BDD-Share-Model" dialog should be hidden
    When user clicks on "BDD-Share-Model" label in gallery
    Then the sharing pane should not list the sharing user

  Scenario: The model is shared with the second account to view and use
    When user picks "Share..." from the context menu of "BDD-Share-Model" label in gallery
    Then "Share BDD-Share-Model" dialog should contain text "Full access"
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDD-Share-Model" dialog
    Then "Send notifications" input in "Share BDD-Share-Model" dialog should be visible
    And "Share BDD-Share-Model" dialog should not contain text "will also be shared"
    When user clicks on OK button in "Share BDD-Share-Model" dialog
    Then the "Share BDD-Share-Model" dialog should close
    And no error or warning balloon should have been shown
    When user clicks on "BDD-Share-Model" label in gallery
    Then the sharing pane should list the sharing user

  Scenario: The owner takes the share back
    When user picks "Share..." from the context menu of "BDD-Share-Model" label in gallery
    Then "Share BDD-Share-Model" dialog should contain text "Full access"
    When user removes the sharing user from "Share BDD-Share-Model" dialog
    And user clicks on OK button in "Share BDD-Share-Model" dialog
    Then the "Share BDD-Share-Model" dialog should close
    And no error or warning balloon should have been shown
    When user clicks on "BDD-Share-Model" label in gallery
    Then the sharing pane should not list the sharing user
    And no errors should have been logged
