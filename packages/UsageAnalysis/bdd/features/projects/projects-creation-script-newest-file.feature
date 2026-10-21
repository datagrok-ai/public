@journey @realizes:views.projects
Feature: A project whose table comes from the user's own script picks up the newest file on reopen
  A JavaScript script reads a folder in My files and returns the CSV file with the highest number at
  the end of its name. Its result is saved as a project with Data sync, the script being the table's
  creation script; a file with a higher number is added to the folder, and the reopened project shows
  the new file's rows, not the rows it was saved with. Translated from the TestTrack case
  Projects/custom-creation-scripts-ui.

  The folder, its two files (cars.csv as data1.csv, 30 rows; iris.csv as data2.csv, 150 rows) and the
  script are the scene, not the subject (the md's automation notes say so), so they go in through the
  API; the second file only after the project is saved. The script is run from its address
  (/func/Admin.<name>?run=true) rather than from the Scripts gallery: what the gallery offers is
  scripts-browser's, and its fuzzy search lists the scripts other features make meanwhile.

  The folder, the files, the script and the project are named with the run's time and removed when the
  feature starts and ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDCcsProject{time}" is on the server
    And a copy of the "System:DemoFiles/cars.csv" file is in the user's files as "BDDCcs{time}/data1.csv"
    And a script "BDDCcsScript{time}" is on the server:
      """
      //language: javascript
      //output: dataframe result
      const folder = grok.shell.user.project.name + ':Home/BDDCcs{time}/';
      const csvFiles = await grok.dapi.files.list(folder, false, 'csv');
      if (csvFiles.length === 0)
        throw new Error('No CSV files found in ' + folder);
      const suffix = (name) => { const m = name.match(/(\d+)(?=\.csv$)/); return m ? parseInt(m[1], 10) : -1; };
      csvFiles.sort((a, b) => suffix(a.fileName) - suffix(b.fileName));
      result = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(csvFiles[csvFiles.length - 1].fullPath));
      """

  Scenario: The script's result is saved with a bar chart, the script as its creation script
    When user opens the address "/func/Admin.BDDCcsScript{time}?run=true"
    Then the current view should be a TableView view
    And the table should have 30 rows
    Given the toolbox pane is shown
    When user clicks on "bar chart" icon in toolbox
    Then bar chart viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "BDDCcsScript{time}" project table in "Save project" dialog should be checked
    When user clicks on "Creation script" button in "BDDCcsScript{time}" project table in "Save project" dialog
    Then "BDDCcsScript{time}" project table in "Save project" dialog should contain text ":BDDCcsScript{time}()"
    When user enters "BDDCcsProject{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCcsProject{time}" uploaded' should have been shown
    And "Share BDDCcsProject{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDCcsProject{time}" dialog
    Then the "Share BDDCcsProject{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    And no table should be left in the workspace

  Scenario: With a newer file in the folder, the reopened project shows its rows
    Given a copy of the "System:DemoFiles/iris.csv" file is in the user's files as "BDDCcs{time}/data2.csv"
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCcsProject{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCcsProject{time} gallery card and the project opens
    Then the current view should be a TableView view
    And the table should have 150 rows
    And the table should have been reloaded by data sync
    And bar chart viewer should be visible
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown
