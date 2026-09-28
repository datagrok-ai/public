@journey @serial @realizes:views.projects
Feature: A project whose table comes from a script picks up the newest file on reopen
  A JavaScript script lists the CSV files of a folder in Browse > Files > My files, takes the one
  with the highest numeric suffix in its name and returns it as a table. The script is run from the
  Browse tree (Run...), a viewer is added to its result, and the result is saved as a project with
  Data sync on, so the script becomes the table's creation script. A file with a higher suffix is
  then added to the folder; the project reopened from its Dashboards card must show that file's
  rows, not the rows it was saved with. Translated from the TestTrack case
  Projects/custom-creation-scripts-ui.

  Fixtures substituted: the md mutates the shared System:DemoFiles/chem folder, which is why it was
  kept manual (a destructive change another run would see, and write access the test account
  lacks). Here the feature owns its folder: BDDCreate<time> under the account's own files (My
  files), holding data_1.csv (3 rows) and later data_2.csv (5 rows), which tell apart by their
  "source" column. The script reads that folder instead of System:DemoFiles/chem and returns the
  file's rows, as the md's script does; its code is saved through the JS API (the script editor is
  not the subject), the run and everything after it go through the UI.

  The folder with its files, the script and the project (with its table and view) are removed when
  the feature starts and when it ends. @serial: the save uploads a table and the reopen reruns the
  script on the server's files.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDCreate{time}" is on the server
    And no folder "BDDCreate{time}" is in the user's files
    And the file "data_1.csv" of the folder "BDDCreate{time}" in the user's files holds:
      """
      source,value
      data_1,1
      data_1,2
      data_1,3
      """
    And a script "BDDCreate{time}" is on the server:
      """
      //language: javascript
      //output: dataframe result
      const folder = grok.shell.user.project.name + ':Home/BDDCreate{time}/';
      const csvFiles = await grok.dapi.files.list(folder, false, 'csv');
      if (csvFiles.length === 0)
        throw new Error('No CSV files found in ' + folder);
      const suffix = (name) => { const m = name.match(/(\d+)(?=\.csv$)/); return m ? parseInt(m[1], 10) : -1; };
      csvFiles.sort((a, b) => suffix(a.fileName) - suffix(b.fileName));
      result = DG.DataFrame.fromCsv(await grok.dapi.files.readAsText(csvFiles[csvFiles.length - 1].fullPath));
      """

  Scenario: The script run from Browse returns the highest-suffixed file
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    Then the "Scripts" view should be current
    When user enters "BDDCreate{time}" into gallery search
    And user picks "Run..." from the context menu of BDDCreate{time} gallery card
    Then the "result" view should be current
    And table "result" should have 3 rows
    And every value of "source" column should match "^data_1$"
    # the Scripts search is the account's setting: emptied again for whoever opens the view next
    Given user opens the Scripts view
    And user switches to the "result" table view
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The result is saved with Data sync, the script as its creation script
    When user clicks on "bar chart" icon in toolbox
    Then the open tableview should have 1 bar chart viewer
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "result" project table in "Save project" dialog should be switched on
    When user clicks on "Creation script" button in "result" project table in "Save project" dialog
    Then "result" project table in "Save project" dialog should contain text "BDDCreate{time}"
    When user enters "BDDCreate{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCreate{time}" uploaded' should have been shown
    And 1 project named "BDDCreate{time}" should be on the server
    And the "result" table of the "BDDCreate{time}" project should be saved with data sync
    And "Share BDDCreate{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDCreate{time}" dialog
    Then the "Share BDDCreate{time}" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors but the project preview's should have been logged

  Scenario: A newer file in the folder is what the reopened project shows
    Given the file "data_2.csv" of the folder "BDDCreate{time}" in the user's files holds:
      """
      source,value
      data_2,10
      data_2,20
      data_2,30
      data_2,40
      data_2,50
      """
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCreate{time}" into gallery search
    And user double-clicks on BDDCreate{time} gallery card
    Then table "result" should have been reloaded by data sync with 5 rows
    And every value of "source" column should match "^data_2$"
    And the open tableview should have 1 bar chart viewer
    And no errors should have been logged
    And no error or warning balloon should have been shown
