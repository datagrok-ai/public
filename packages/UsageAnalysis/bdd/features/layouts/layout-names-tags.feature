@realizes:views.layouts @journey
Feature: A layout saved from the toolbox gets a name and tags
  Save in the Layouts pane of the toolbox saves the layout of the table view and opens it in the
  context panel, named after the table. The name typed into the panel renames it everywhere the UI
  shows it, the pencil next to Tags adds a tag, the filter of the pane finds the layout by part of
  its name and by the tag, and once the table is opened again the card shows what the server keeps
  and a click on it applies the layout.

  A click on the Layouts header blocks every change of the current object for 2 s, which would keep
  the saved layout out of the panel; the click on a column header after it releases that block and
  makes the column current, which the panel shows before Save. Collapsing and expanding the section
  does not fetch the cards again, so the server's copy is read by closing the table and opening it
  anew. The filter applies 200 ms after typing: each search that must keep the card follows one
  that hides it. A card keeps the copy of the layout it was built with: Edit tags... is tried on a
  second layout, saved and renamed after the table was opened again, so the first one's rename
  is read back from the server before that.

  Background:
    Given user is logged in
    And the layouts named "bdd-layout-{time}" are deleted when the feature ends
    And the layouts named "BDD height vs weight-{time}" are deleted when the feature ends
    And the layouts named "BDD tags-{time}" are deleted when the feature ends
    And user opens demog dataset keeping the first 100 rows as "bdd-layout-{time}"
    And user adds a scatter plot viewer with:
      | xColumnName | HEIGHT |
      | yColumnName | WEIGHT |
    And the toolbox pane is shown
    And the context panel is open
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"

  Scenario: Save in the Layouts pane names the layout after the table and opens it in the context panel
    When user clicks on Save button in layouts pane
    Then "bdd-layout-{time}" layout card should be visible
    And the context panel should show "bdd-layout-{time}"
    And the title of context panel should be "bdd-layout-{time}"
    And context panel should contain text "Layout saved"
    And no errors should have been logged

  Scenario: The name typed into the context panel renames the layout on its card and in the panel
    When user types "BDD height vs weight-{time}" into Name details field in context panel
    And user presses Enter
    Then "BDD height vs weight-{time}" layout card should be visible
    And "bdd-layout-{time}" layout card should be absent
    And value of Name details field in context panel should have text "BDD height vs weight-{time}"
    And the context panel should show "BDD height vs weight-{time}"
    And no errors should have been logged

  @known-failure
  Scenario: The context panel title follows the new name (GROK-21137)
    Then the title of context panel should be "BDD height vs weight-{time}"

  Scenario: The pencil next to Tags adds a tag that the card shows
    Given the context panel should show "BDD height vs weight-{time}"
    When user hovers over Tags details row in context panel
    And user clicks on edit icon of Tags details row in context panel
    And user types "anthropometry" at the caret
    And user presses Enter
    Then "BDD height vs weight-{time}" layout card should contain text "#anthropometry"
    And no errors should have been logged

  Scenario: The filter of the Layouts pane finds the layout by part of its name and by the tag
    When user types "bdd-layout-{time}" into layouts filter
    Then there should be 0 visible layout cards
    When user types "height vs weight-{time}" into layouts filter
    Then "BDD height vs weight-{time}" layout card should be visible
    And there should be 1 visible layout card
    When user types "bdd-layout-{time}" into layouts filter
    Then there should be 0 visible layout cards
    When user types "#anthropometry" into layouts filter
    Then "BDD height vs weight-{time}" layout card should be visible
    When user clears layouts filter
    Then "BDD height vs weight-{time}" layout card should be visible

  Scenario: The table opened again lists the renamed, tagged layout, and a click on its card applies it
    When user closes all views
    And user opens demog dataset keeping the first 100 rows as "bdd-layout-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    Then "BDD height vs weight-{time}" layout card should be visible
    And "bdd-layout-{time}" layout card should be absent
    And "BDD height vs weight-{time}" layout card should contain text "#anthropometry"
    And the open tableview should have 0 scatter plot viewers
    When user clicks on "BDD height vs weight-{time}" layout card
    Then the open tableview should have 1 scatter plot viewer
    And "xColumnName" property of scatter plot viewer should be "HEIGHT"
    And no errors should have been logged

  Scenario: A second layout saved and renamed takes a tag through Edit tags... on its card
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    When user clicks on Save button in layouts pane
    Then the context panel should show "bdd-layout-{time}"
    When user types "BDD tags-{time}" into Name details field in context panel
    And user presses Enter
    Then "BDD tags-{time}" layout card should be visible
    When user picks "Edit tags..." from the context menu of "BDD tags-{time}" layout card
    Then "Edit tags" dialog should be visible
    When user types "anthropometry2" into Tags input in "Edit tags" dialog
    And user clicks on OK button in "Edit tags" dialog
    Then "Edit tags" dialog should be absent

  @known-failure
  Scenario: Edit tags... keeps the name typed into the context panel (GROK-21138)
    Then "BDD tags-{time}" layout card should be visible
    And "bdd-layout-{time}" layout card should be absent
