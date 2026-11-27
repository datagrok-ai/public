@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). demog.csv is opened off camera from Browse > Files > Demo by
  a double click, so its table carries a creation script, and a scatter plot is added. On camera the ribbon's Save
  opens the Save project dialog, whose row for the table shows Data sync on: the switch's tooltip says what it does,
  and Creation script unfolds the script that reads the file again. Data sync is switched off and the dashboard is
  saved as "Demographics snapshot", its table stored as a snapshot; the Share offer the first save opens is
  cancelled.

  The project, with its table, view and picture, is deleted when the feature starts and ends.

  Scenario: Where a dashboard gets its data
    Data sync decides whether a dashboard reads its data again from the source or keeps a saved copy.

    Given user is logged in
    And simple mode is off
    And no project named "Demographics snapshot" is on the server
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And user double-clicks Files---Demo---demog.csv tree node inside browse tree
    And the "demog" view should be current
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE    |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "This table came from a file. A dashboard can read it again or keep a copy"
    # caption: SAVE lists the table with its Data sync switch
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    # caption: Data sync on: each time it opens, the dashboard reads the file again
    When user hovers over "Data sync" text in "demog" project table in "Save project" dialog
    # caption: Creation script: the steps that make the table again
    And user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    # caption: Data sync off: the data is saved with the dashboard, as it is now
    And user switches off Data sync switch in "demog" project table in "Save project" dialog
    Then Data sync switch in "demog" project table in "Save project" dialog should be switched off
    When user enters "Demographics snapshot" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the "demog" table of the "Demographics snapshot" project should be saved as a snapshot
    And "Share Demographics snapshot" dialog should be visible
    When user clicks on CANCEL button in "Share Demographics snapshot" dialog
    Then the "Share Demographics snapshot" dialog should close
    When the video says "Keep it on for current data, off for a report that should not change"
    And user hovers over scatter plot viewer
    Then no errors should have been logged
