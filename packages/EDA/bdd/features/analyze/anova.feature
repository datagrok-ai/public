@eda @realizes:ml.menu.analyze.group-comparison.anova
Feature: One-way ANOVA
  ML | Analyze | Group Comparison | ANOVA... over demog with the dialog's own choices: AGE across
  RACE at alpha 0.05. Translated from files/TestTrack/EDA/anova.md and the package's
  playwright/anova.test.ts.

  The case still names an Analysis and an F-test tab; the platform replaced them (GROK-20194) with
  the conclusion on the box plot's description and a docked table of the test, which is what this
  feature reads. The p-value is the one demog gives, so a change in the computation fails here.

  Background:
    Given user is logged in
    And user opens demog dataset

  Scenario: The dialog opens on a category, a feature and a significance level
    When user picks "ML > Analyze > Group Comparison > ANOVA..." from the top menu
    Then "ANOVA" dialog should be visible
    And editor of Category input in "ANOVA" dialog should have text "RACE"
    And editor of Feature input in "ANOVA" dialog should have text "AGE"
    And Alpha input in "ANOVA" dialog should have value "0.05"
    And Run button in "ANOVA" dialog should be enabled

  Scenario: Running it docks a box plot with the conclusion and the table of the test
    When user picks "ML > Analyze > Group Comparison > ANOVA..." from the top menu
    And user clicks on Run button in "ANOVA" dialog
    Then the top menu command should have completed
    And "ANOVA" dialog should be hidden
    And box plot viewer should be visible
    And "Description" property of box plot viewer should contain "doesn't affect"
    And "Description" property of box plot viewer should contain "p = 0.176"
    And box plot viewer should be painted
    And table "ANOVA result" should be open
    And table "ANOVA result" should have 1 row
    And table "ANOVA result" should have columns "Conclusion, Source of variance, F, df₁, df₂, F-critical, p-value"
    And second grid viewer should be bound to table "ANOVA result"
    And no error or warning balloon should have been shown
    And no errors should have been logged
