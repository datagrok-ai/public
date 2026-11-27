@guide @tutorial
Feature: How-to · Scripting
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The console, opened from the status bar's Console toggle, logs
  what the user does as commands: the Add New Column dialog (BMI from WEIGHT and HEIGHT) comes out as an
  AddNewColumn("demog", <formula>, "BMI") line. Up in the console's command line brings the last command back, and
  the same command with another formula and name, entered there, adds a second column (pounds) without the dialog.
  demog is opened off camera with the toolbox, the context panel and Browse hidden. Nothing is made on the server:
  the two columns live only in the open table.

  Scenario: The console
    The console logs what you do as commands. Run a command again from there, with other values.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: The console: commands and their results
    When user clicks on "Console" status bar toggle
    Then console should be visible
    When user notes the console output
    # caption: Add a column with the dialog, as usual
    And user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "BMI" into column name input
    And user replaces the code of code editor in "Add New Column" dialog with "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"
    And user clicks on OK button in "Add New Column" dialog
    Then the table should have a column "BMI"
    # caption: The console logged it as a command
    And the console should show "AddNewColumn(\"demog\", \"${WEIGHT} / Pow(${HEIGHT} / 100, 2)\", \"BMI\")"
    # caption: Type it here with another formula and name, then Enter
    When user clicks on console command line
    And user pastes "AddNewColumn(\"demog\", \"${WEIGHT} * 2.2\", \"Weight in lb\")" into console command line
    And user presses Enter
    Then the table should have a column "Weight in lb"
    # 73.20 kg * 2.2
    And the "Weight in lb" cell of row 1 should be displayed as "161.04"
    # caption: A command's name alone opens a form for its inputs
    When user types "SelectRandomRows" into console command line
    And user presses Enter
    Then "selectRandomRows" dialog should be visible
    When user enters "100" into Count input in "selectRandomRows" dialog
    And user clicks on OK button in "selectRandomRows" dialog
    Then 100 rows should be selected
    # caption: The console shows the command with the values you chose
    And the console should show "selectRandomRows(\"demog\", 100, 42, true)"
    And no errors should have been logged
