@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A file from the computer uploaded to the account's own folder,
  Browse > Files > My files: My files is clicked (its folder view), the file (files/solubility-assay.csv of this
  folder, 16 rows) is dropped on the folder's view, which uploads it into My files; it shows in the folder view and
  in the tree, a click previews it and a double-click opens it as a table. Dropping a file anywhere else opens it
  without uploading (01-getting-started/02). The uploaded file is deleted when the feature ends (and a leftover one
  when it starts).

  Scenario: Upload a file to your files
    Drop a file from your computer into My files: it is kept on the server and opens from Browse.

    Given user is logged in
    And simple mode is off
    And no file "solubility-assay.csv" is in the user's files, now or when the feature ends
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user clicks on Files---My-files tree node inside browse tree
    Then the "My files" view should be current
    When the video says "Browse > Files > My files: your own folder on the server"
    And the video says "Drag a file from your computer onto it"
    And user drops the "features/howto/14-files/files/solubility-assay.csv" file of the project onto gallery to upload it
    Then solubility-assay.csv link in gallery should be visible
    And the file "solubility-assay.csv" should be in the user's files on the server
    When the video says "It is uploaded: in the folder and in the tree"
    And user hovers over solubility-assay.csv link in gallery
    And user clicks on Files---My-files---solubility-assay.csv tree node inside browse tree
    Then the "solubility-assay" view should be current
    And no table should be left in the workspace
    When the video says "One click previews it"
    And user double-clicks on Files---My-files---solubility-assay.csv tree node inside browse tree
    Then the table should have 16 rows
    And no errors should have been logged
    When the video says "A double-click opens it as a table"
