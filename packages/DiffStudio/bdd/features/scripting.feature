@journey @diffstudio @realizes:diffstudio.app.diff-studio
Feature: The equations behind a model, and the script they become
  The Edit toggle opening the equations editor over a model, and the </> command turning the model
  into a script view. Translated from files/TestTrack/DiffStudio/scripting.md and the spec beside it.

  The chart of the script view is the one the manual checklist gave up on (M-1.5): its canvas
  answers toDataURL with a blank image. It is a real line chart of the platform, docked under a tab
  of the function view the script runs in, so the redraw is the viewer's own word — once its tab is
  the one shown, since the viewers of the other tabs stay in the DOM without a rectangle. The table
  is the grid under the other tab, and its rows follow the model's Final time.

  The tail of the case — tag the script as a model, save it, and open it again from the Model Hub —
  is here too (M-1.6 of the manual checklist). Saving leaves a script on the stand, so the step that
  saves notes which scripts were there first and deletes exactly the one it added when the feature
  ends; a model-tagged script that outlived a run would show up in everyone's Model Hub.

  Background:
    Given user is logged in
    And user opens the "Bioreactor" model of the Diff Studio library

  Scenario: The model shows its inputs before anything is edited
    Then the "Bioreactor" view should be current
    And "Process mode" input should be visible
    And code editor should be absent

  Scenario: Edit opens the equations editor in place of the form
    When user clicks on Edit ribbon item
    Then code editor should be visible
    And "Process mode" input should be absent

  Scenario: The angle brackets turn the model into a script
    When user clicks on "</>" ribbon item
    Then Sensitivity ribbon item should be absent
    And "Run script (F5)" icon should be visible

  Scenario: The script is tagged as a model and saved
    When user puts "//tags: model" on the first line of code editor
    Then code editor should contain the text "//tags: model"
    When user saves the script

  Scenario: The script runs, and its table and chart follow the input it exposes
    When user clicks on "Run script (F5)" icon
    Then Final input should be visible
    And "Bioreactor / Grid" tab should be selected
    And the "rows" reading of grid viewer should be "1001"
    When user clicks on "Bioreactor / DiffStudio Facet" tab
    And user takes a snapshot of line chart viewer
    And user enters "500" into Final input
    Then Final input should have value "500"
    And line chart viewer should have repainted
    When user clicks on "Bioreactor / Grid" tab
    Then the "rows" reading of grid viewer should be "501"
    And no errors should have been logged

  Scenario: The Model Hub lists the saved script, and it answers to its inputs there
    Given user opens the Model Hub
    Then the Model Hub should list the saved script
    When user opens the saved script from the Model Hub
    Then Final input should be visible
    And "Bioreactor / Grid" tab should be selected
    When user clicks on "Bioreactor / DiffStudio Facet" tab
    And user takes a snapshot of line chart viewer
    And user enters "800" into Final input
    Then Final input should have value "800"
    And line chart viewer should have repainted
    When user clicks on "Bioreactor / Grid" tab
    Then the "rows" reading of grid viewer should be "801"
    And no errors should have been logged

  Scenario: Refresh re-fetches the catalog, so a model removed behind its back disappears
    When the saved script is deleted on the server
    And user clicks on model hub refresh icon
    Then the Model Hub should not list the saved script
    And no errors should have been logged
