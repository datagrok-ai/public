@regression
Feature: Activity Cliffs on a table with no numeric column
  smiles_only holds 1000 molecules and no other column, so the Activity Cliffs dialog has nothing to
  offer as the activity. OK is disabled from the moment the dialog opens: a click that lands on it
  at once runs nothing, logs no error and leaves the dialog open (GROK-20915 — OK used to be
  clickable for the first seconds, and that click failed with a TypeError balloon).

  The click is the pointer landing on the button as soon as the button is on screen, not the
  library's click, which waits for an enabled target and so would wait out the very seconds the
  defect lived in. OK's state is recorded in the task it enters the page and again as the click
  lands, together with whether the click point is on OK (a dialog still laying itself out moves it),
  and the absence of a balloon and of an error is held over two seconds, not read once.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens smiles-only dataset

  @realizes:GROK-20915
  Scenario: OK clicked at once runs nothing and the dialog stays open
    Given user watches the OK button of the next dialog
    When user picks "Chem > Analyze > Activity Cliffs..." from the top menu
    And user clicks on OK button in "Activity Cliffs" dialog at once
    Then the OK button should have been disabled when it appeared and when it was clicked
    And no error or warning balloon and no error should appear for 2 seconds
    And "Activity Cliffs" dialog should be visible
    And OK button in "Activity Cliffs" dialog should be disabled
    And no new column should have been added
    And the table should have 1000 rows
