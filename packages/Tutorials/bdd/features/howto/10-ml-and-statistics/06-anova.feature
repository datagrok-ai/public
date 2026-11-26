@guide @tutorial
Feature: How-to · ML
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). ML > Analyze > Group Comparison > ANOVA... is EDA's one-way
  ANOVA (Welch by default), computed in the browser. On demog-1000 the dialog starts on RACE and AGE; the clip picks
  DIS_POP as the category (six disease groups, mean ages 37 to 52: F about 43 in a local check), keeps AGE and alpha
  0.05, and runs it: a box plot of AGE by DIS_POP docks on the right with the conclusion on top ("DIS_POP" affects
  the "AGE", p < 0.001), and the "ANOVA result" table docks under it (Conclusion, F, degrees of freedom, F-critical,
  p-value). demog-1000 is opened off camera; nothing is made on the server.

  Scenario: Compare groups with ANOVA
    One-way ANOVA tells whether a number differs between groups, and the box plot shows how.

    Given user is logged in
    And simple mode is off
    And user opens demog-1000 dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "ANOVA checks whether age differs between six disease groups" at the bottom
    # caption: ANOVA is under ML > Analyze > Group Comparison
    And user picks "ML > Analyze > Group Comparison > ANOVA..." from the top menu
    Then "ANOVA" dialog should be visible
    # caption: Category: the column that splits rows into groups
    When user selects "DIS_POP" in Category input in "ANOVA" dialog
    Then editor of Category input in "ANOVA" dialog should have text "DIS_POP"
    # caption: Feature: the number to compare between the groups
    When user hovers over editor of Feature input in "ANOVA" dialog
    Then editor of Feature input in "ANOVA" dialog should have text "AGE"
    # caption: Alpha: how sure you want to be, 0.05 is usual
    When user hovers over Alpha input in "ANOVA" dialog
    Then Alpha input in "ANOVA" dialog should have value "0.05"
    When user clicks on Run button in "ANOVA" dialog
    Then "ANOVA" dialog should be hidden
    And box plot viewer should be visible
    And description of box plot viewer should contain text "affects"
    And table "ANOVA result" should have 1 row
    When the video says "The answer is on top of the box plot, with its p-value" at the bottom
    And user hovers over description of box plot viewer
    And the video says "Below: F above F-critical, p below alpha, so the groups differ" at the top
    And user hovers over second grid viewer
    Then second grid viewer should be bound to table "ANOVA result"
    And no error or warning balloon should have been shown
    And no errors should have been logged
