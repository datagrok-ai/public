@guide
Feature: Open a JMP file in the browser
  A guide: the answer to "how do I open a .jmp file in Datagrok, without JMP?". Click it in the
  Browse tree: the table opens with the analyses saved in it as viewers, and the "About this JMP
  file" icon tells what the file is and what didn't come across. Demo: the published lifecycle DoE
  of Kunzelmann et al. 2024 (Bioengineering 11:1089, CC BY 4.0), 21 runs saved by JMP 18.1.0; its
  "Fit OLS" analysis opens as three Actual by Predicted plots, with the file's spec limits as lines.

  Scenario: Double-click a JMP file to open it with its saved analyses
    Given user is logged in
    And simple mode is off
    And no table view is open
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---JmpReader tree node inside browse tree is expanded
    And Files---App-Data---JmpReader---public tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---JmpReader---public---mdpi-bioengineering-11-01089-wp1.jmp tree node inside browse tree
    Then the "mdpi-bioengineering-11-01089-wp1" table view should open with 21 rows
    And the open tableview should have 3 scatter plot viewers
    When user clicks on "About this JMP file" icon
    Then "Saved analyses" section should be visible
    And "Saved analyses" section should contain text "Scatter plot"
    When user clicks on CANCEL button
    Then the "mdpi-bioengineering-11-01089-wp1.jmp" dialog should close
