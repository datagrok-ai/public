@journey @realizes:views.functions @realizes:views.projects
Feature: A function run straight from its URL
  /func/<Namespace>.<Name>?<param>=<value> opens a script's parameter form; &run=true runs it at
  once and shows only the result; the copy icon left of RUN gives such a link for the current
  values; a run link without a required value falls back to the form and names what is missing.
  Two JavaScript scripts (demog(rows) with a default, and one with a required rows) are made
  through the API (writing a script is the scripts features' subject), and the namespace of the address is read from the Grok name the first script's
  Links... dialog shows. Translated from the TestTrack case Projects/function-url-run.

  The addresses are opened in the same tab and typed (the md pastes them into a new one), each after
  a Close All, as a new tab starts with nothing open: an address of a function whose view is open
  does not open it again. Kept without (see the request document): opening the run link from the
  clipboard — the clipboard is claimed to hold it, and the same address is typed.

  A run link without the required value (/func/Admin.<script>?label=%22y%22&run=true — a string
  value in the address is quoted, as the copy icon writes it) falls back to the form with an
  "Unable to run" warning naming rows (core#446). The warning is claimed in the scenario that opens
  the address, and once: a scenario starts with the balloons cleared, and a matched balloon is taken.

  The scripts are named with the run's time (letters and digits only) and removed when the feature
  starts and ends.

  Background:
    Given user is logged in
    And the context panel is open
    And the browse panel is open
    And a script "BDDUrlRunScript{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 10
      //output: dataframe df
      df = grok.data.demo.demog(rows);
      """
    And a script "BDDUrlRunRequired{time}" is on the server:
      """
      //language: javascript
      //input: int rows
      //input: string label = "x"
      //output: dataframe df
      df = grok.data.demo.demog(rows);
      """

  Scenario: The script's Links... dialog gives its Grok name
    When user opens the address "/scripts"
    Then the "Scripts" view should be current
    When user clears gallery search
    And user types "BDDUrlRunScript{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on "BDDUrlRunScript{time}" link in gallery
    Then the context panel should show "BDDUrlRunScript{time}"
    When user clicks on "Links..." link in Details pane in context panel
    Then "Links to BDDUrlRunScript{time}" dialog should be visible
    And "Grok name" input in "Links to BDDUrlRunScript{time}" dialog should have value "Admin:BDDUrlRunScript{time}"
    When user presses Escape
    Then the "Links to BDDUrlRunScript{time}" dialog should close
    # the gallery keeps its search text for the next visit: leave it empty
    When user clears gallery search

  Scenario: Without run the form opens with the URL's value
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/func/Admin.BDDUrlRunScript{time}?rows=25"
    Then "Rows" input should have value "25"
    And RUN button should be visible
    And status bar should contain text "Rows: 0"

  Scenario: The copy icon gives a run link for the current value
    When user enters "40" into "Rows" input
    And user hovers over copy icon
    Then tooltip should contain text "Copy a link that runs this function with the current parameters"
    When user clicks on copy icon
    Then the clipboard should contain text "/func/Admin.BDDUrlRunScript{time}?rows=40&run=true"

  Scenario: The run link runs at once and shows only the result
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/func/Admin.BDDUrlRunScript{time}?rows=40&run=true"
    Then the current view should be a TableView view
    And the table should have 40 rows
    And status bar should contain text "Rows: 40"
    And the page address should contain "run=true"
    And RUN button should be absent
    Given the toolbox pane is shown
    Then "Rows" input in Source pane in toolbox should have value "40"
    And REFRESH button in Source pane in toolbox should be visible

  Scenario: A value changed in the result view is refreshed and carried by the link
    When user enters "15" into "Rows" input in Source pane in toolbox
    And user clicks on REFRESH button in Source pane in toolbox
    Then the "rows" reading of grid should be 15
    And status bar should contain text "Rows: 15"
    When user hovers over copy icon in Source pane in toolbox
    Then tooltip should contain text "?rows=15&run=true"

  Scenario: run=false keeps the form
    When user moves the pointer away from Source pane in toolbox
    And user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/func/Admin.BDDUrlRunScript{time}?rows=25&run=false"
    Then "Rows" input should have value "25"
    And RUN button should be visible
    And status bar should contain text "Rows: 0"

  @realizes:core#446
  Scenario: The run link without a required value falls back to the form and names it (core#446)
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/func/Admin.BDDUrlRunRequired{time}?label=%22y%22&run=true"
    Then the "BDDUrlRunRequired{time}" view should be current
    And a warning balloon containing 'Unable to run "BDDUrlRunRequired{time}": rows' should have been shown
    And "Label" input should have value "y"
    And RUN button should be visible
    And status bar should contain text "Rows: 0"
    When user enters "5" into "Rows" input
    And user clicks on RUN button
    Then the current view should be a TableView view
    And the table should have 5 rows
