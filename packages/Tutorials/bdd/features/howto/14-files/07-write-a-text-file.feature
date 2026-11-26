@guide @tutorial
Feature: How-to · Files
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A text file written in Browse > Files > My files: Create file...
  from the folder view's context menu, the name typed, the new file clicked (a text file previews in an editor,
  `file_editors.dart` _txtFileEditor), a few lines typed, and the ribbon's SAVE, which shows once the text has
  changed, writes it to the server. The file is deleted when the feature ends (and a leftover one when it starts).

  Scenario: Write a text file
    Keep notes next to your data: make a text file in your files, type, and save.

    Given user is logged in
    And simple mode is off
    And no file "assay-notes.txt" is in the user's files, now or when the feature ends
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user clicks on Files---My-files tree node inside browse tree
    Then the "My files" view should be current
    When the video says "Right-click the folder view to make a new file"
    And user picks "Create file..." from the context menu of gallery
    Then "Create file" dialog should be visible
    When user enters "assay-notes.txt" into Name input in "Create file" dialog
    And user clicks on OK button in "Create file" dialog
    Then the "Create file" dialog should close
    And assay-notes.txt link in gallery should be visible
    When user expands Files---My-files tree node inside browse tree
    And user clicks on Files---My-files---assay-notes.txt tree node inside browse tree
    Then the "assay-notes.txt" view should be current
    When the video says "A text file opens in an editor"
    And user clicks on code editor
    And user types "Solubility assay, plate 3: re-run CMP-106 at pH 6.5" at the caret
    Then SAVE button should be visible
    When the video says "SAVE shows once the text changes"
    And user clicks on SAVE button
    Then the file "assay-notes.txt" in the user's files should hold the text "re-run CMP-106 at pH 6.5"
    And SAVE button should be hidden
    And no errors should have been logged
    When the video says "Saved to your files on the server"
