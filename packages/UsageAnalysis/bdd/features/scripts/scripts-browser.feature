@journey @full-stand @serial @realizes:views.scripts
Feature: The Scripts view and a script's context panel
  The Scripts view finds a script by search; clicked, the script fills the context panel with its
  details (inputs and outputs), its activity (the creation and every run), its sharing and its chats;
  Share... gives the second account access, a chat message stays on the script; the gallery switches
  its view mode and its order; a shared sample script (ACF) runs on TSLA and opens in the editor.
  Translated from files/TestTrack/Scripts/browser.md and playwright-public/scripts/scripts-browser.test.ts.

  The script is this feature's own ({time} in its name); it is deleted with its chat at the end,
  and both are checked gone. R (the script and ACF) runs in a container, so the
  feature is @full-stand.

  Not translated, and why: the md's "Usage" tab is the Activity pane now (History stays at 0 for a
  run started this way), and its entries are not claimed at all: the pane
  expands under a hand click and loads "<user> ran <script>", but in this journey it stays collapsed
  whatever the click lands on, and a claim on a collapsed pane would be a claim on its header. The
  count after a run is read in a fresh Scripts view: the context panel renders the current object
  once, and the same card clicked again keeps the count it showed. A run through the JS API is not
  counted there at all, and the platform records a script made through it with a delay of its own,
  so the claim is "at least one" — the md's run is the UI's. The Sort list menu is not claimed here: after a
  view-mode switch the gallery reloads for longer than a click waits, and the Users suite already
  holds that menu and the order it gives. The old spec checked "Activity runs after >= before",
  which an unchanged count passed.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And a script "BddScriptBrowser{time}" is on the server:
      """
      #language: r
      #input: dataframe table
      #output: int count
      #output: string newParam
      count <- nrow(table) * ncol(table)
      newParam <- "test"
      """
    And the context panel is open
    And user opens the Scripts view

  Scenario: The view modes and the order of the gallery
    When user clicks on "Switch to brief view" icon inside gallery toolbar
    Then the gallery should be in brief mode
    When user clicks on "Switch to grid view" icon inside gallery toolbar
    Then the gallery should be in grid mode
    When user clicks on "Switch to card view" icon inside gallery toolbar
    Then the gallery should be in card mode
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A search narrows the gallery to the script
    When user types "BddScriptBrowser{time}" into gallery search
    Then gallery counter should have text "1"
    And "BddScriptBrowser{time}" link in gallery should be visible
    And no errors should have been logged

  Scenario: The script fills the context panel
    When user clicks on "BddScriptBrowser{time}" link in gallery
    Then the context panel should show "BddScriptBrowser{time}"
    And the following elements should be visible:
      | "Details" section in context panel  |
      | "Script" section in context panel   |
      | "Run" section in context panel      |
      | "Sharing" section in context panel  |
      | "Chats" section in context panel    |
      | "Activity" section in context panel |
    When user clicks on "Details" pane header in context panel
    Then "Details" section in context panel should contain text "count, newParam"
    And "Details" section in context panel should contain text "table"
    And no errors should have been logged

  Scenario: Activity records what the script has been through
    Given user opens cars dataset
    And user opens the Scripts view
    When user types "BddScriptBrowser{time}" into gallery search
    And user picks "Run..." from the context menu of "BddScriptBrowser{time}" link in gallery
    And user selects "cars" in Table input in "BddScriptBrowser{time}" dialog
    And user clicks on OK button in "BddScriptBrowser{time}" dialog
    Then the "BddScriptBrowser{time}" dialog should close
    # the panel renders the current object once: a fresh Scripts view shows the script anew
    When user closes the current view
    And user opens the Scripts view
    And user types "BddScriptBrowser{time}" into gallery search
    And user clicks on "BddScriptBrowser{time}" link in gallery
    Then the context panel should show "BddScriptBrowser{time}"
    And the "Activity" pane of the context panel should count at least 1
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Share... gives the second account access
    Then the sharing pane should not list the sharing user
    When user picks "Share..." from the context menu of "BddScriptBrowser{time}" link in gallery
    Then "Share BddScriptBrowser{time}" dialog should be visible
    # the dialog fetches the entity's project after it opens; OK before that throws "Not initialized"
    And "Share BddScriptBrowser{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BddScriptBrowser{time}" dialog
    And user clicks on OK button in "Share BddScriptBrowser{time}" dialog
    Then the "Share BddScriptBrowser{time}" dialog should close
    When user clicks on "BddScriptBrowser{time}" link in gallery
    Then the sharing pane should list the sharing user
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A chat message stays on the script
    When user clicks on "Chats" pane header in context panel
    When user types "bdd chat {time}" into chat input in context panel
    And user presses Enter in chat input in context panel
    Then "Chats" section in context panel should contain text "bdd chat {time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The ACF sample runs on TSLA's Close
    Given the "System:DemoFiles/TSLA.csv" file is loaded as a table
    And user opens the Scripts view
    When user types "ACF" into gallery search
    And user notes the console output
    And user picks "Run..." from the context menu of "ACF" link in gallery
    Then "ACF" dialog should be visible
    When user selects "TSLA" in Data input in "ACF" dialog
    And user clicks on editor of Columns input in "ACF" dialog
    Then "Select columns..." dialog should be visible
    And the "text of cell 5 of __name" reading of grid viewer in "Select columns..." dialog should be "Close"
    When user clicks on the "cell 5 of x" area of grid viewer in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "1 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Columns input in "ACF" dialog should contain text "(1)"
    When user clicks on OK button in "ACF" dialog
    Then the "ACF" dialog should close
    And the console should show "ACF("
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Edit... opens the ACF sample in the editor
    Given user opens the Scripts view
    When user types "ACF" into gallery search
    And user picks "Edit..." from the context menu of "ACF" link in gallery
    Then the "ACF" view should be current
    And code editor should contain the text "#name: ACF"
    And no errors should have been logged
    When user closes the current view
    And user opens the Scripts view
