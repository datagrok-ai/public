@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Files kept in order in Browse > Files > My files: a folder made
  from the folder view's context menu (Create folder...), a file dragged onto it in the tree (a move inside one
  share, `files_view.dart` bind), and the file renamed from its own context menu (Rename...). The file
  (files/solubility-assay.csv of this folder) is put into My files off camera; the folder "Solubility" and all it
  holds are deleted when the feature ends (and a leftover one when it starts), the file at its first place too.

  Scenario: Organize your files
    Make folders in My files, drag files into them, and rename them.

    Given user is logged in
    And simple mode is off
    And no file "Solubility" is in the user's files, now or when the feature ends
    And the "features/howto/14-files/files/solubility-assay.csv" file of the project is in the user's files as "solubility-assay.csv"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user clicks on Files---My-files tree node inside browse tree
    Then the "My files" view should be current
    And solubility-assay.csv link in gallery should be visible
    When the video says "Browse > Files > My files: your own folder on the server"
    And the video says "Right-click the folder view to make a folder"
    And user picks "Create folder..." from the context menu of gallery
    Then "Create folder" dialog should be visible
    When user enters "Solubility" into Name input in "Create folder" dialog
    And user clicks on OK button in "Create folder" dialog
    Then the "Create folder" dialog should close
    And Solubility link in gallery should be visible
    When the video says "Drag a file onto a folder to move it there"
    And user expands Files---My-files tree node inside browse tree
    And user drags Files---My-files---solubility-assay.csv tree node inside browse tree to Files---My-files---Solubility tree node inside browse tree
    Then the file "Solubility/solubility-assay.csv" should be in the user's files on the server
    And the file "solubility-assay.csv" should not be in the user's files on the server
    When user expands Files---My-files---Solubility tree node inside browse tree
    Then Files---My-files---Solubility---solubility-assay.csv tree node inside browse tree should be visible
    When the video says "Right-click a file to rename it"
    And user picks "Rename..." from the context menu of Files---My-files---Solubility---solubility-assay.csv tree node inside browse tree
    Then "File name" input in Rename dialog should have value "solubility-assay.csv"
    When user enters "solubility-ph-study.csv" into "File name" input in Rename dialog
    And user clicks on OK button in Rename dialog
    Then the "Rename" dialog should close
    And the file "Solubility/solubility-ph-study.csv" should be in the user's files on the server
    And the file "Solubility/solubility-assay.csv" should not be in the user's files on the server
    And Files---My-files---Solubility---solubility-ph-study.csv tree node inside browse tree should be visible
    And no errors should have been logged
    When the video says "The file has its new name, in its new folder"
