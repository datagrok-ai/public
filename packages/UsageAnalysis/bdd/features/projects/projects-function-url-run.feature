@journey @serial @realizes:views.functions @realizes:views.projects
Feature: A function run straight from its URL
  /func/<Namespace>.<Name>?<param>=<value> opens a script's parameter form; &run=true runs it at
  once and shows only the result; the copy icon left of RUN gives such a link for the current
  values. A JavaScript script, demog(rows) with a default, is written in the Scripts view, and the
  namespace of the address is read from the Grok name the script's Links... dialog shows.
  Translated from the TestTrack case Projects/function-url-run.

  The addresses are opened in the same tab and typed (the md pastes them into a new one), each after
  a Close All, as a new tab starts with nothing open: an address of a function whose view is open
  does not open it again. Kept without (see the request document): opening the run link from the
  clipboard — the clipboard is claimed to hold it, and the same address is typed.

  Not translated: the md's last step, a run link without a required value (its script, with a
  required rows, is not made either). On localhost (1.28.0, f98337ae1b)
  /func/Admin.<script>?label=%22y%22&run=true runs the script with rows empty (a table of 100,000
  rows) instead of falling back to the form with an "Unable to run" warning — a suspected defect
  (GROK-20511 introduced the run links), described outside the repository; the md's unquoted
  label=y is read as an expression ("Variable "y" not found").

  The script is named with the run's time (letters and digits only) and removed when the feature
  starts and ends. It is serial: the Scripts view's search is shared with the features that use it.

  Background:
    Given user is logged in
    And the browse panel is open
    And no script named "BDDUrlRunScript{time}" is on the server

  Scenario: The script is written in the Scripts view
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    Then the "Scripts" view should be current
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    When user replaces the code of code editor with "//name: BDDUrlRunScript{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//input: int rows = 10" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "df = grok.data.demo.demog(rows);" to code editor
    And user saves the script
    Then 1 script named "BDDUrlRunScript{time}" should be on the server

  Scenario: The script's Links... dialog gives its Grok name
    When user closes the current view
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
