@journey @viewers @realizes:viewers.box-plot
Feature: Box plot group comparison
  The group-comparison ladder: the bare p-value and its reveal icon, the overall test and the
  on-chart method selector, one-way ANOVA on three groups, a control group and its comparisons
  table, two-way ANOVA, closing the comparison; then covariate adjustment — regress-out, ratio,
  ANCOVA with a control — and the matched per-stratum baseline with its Simpson's-paradox cue.
  One journey on demog-1000 with a box plot of AGE by SEX; every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE |
      | Category 1 | SEX |

  Scenario: The bare p-value and its reveal icon
    Then "Show P Value" property of box plot viewer should be "true"
    And "Show Group Comparison" property of box plot viewer should be "false"
    And box plot viewer should have a "p value" area
    When user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "t-test"
    And show-group-stats icon in box plot viewer should be visible
    When user clicks on show-group-stats icon in box plot viewer
    Then "Show Group Comparison" property of box plot viewer should be "true"
    And box plot viewer should have a "group comparison" area

  Scenario: The overall test and the method selector
    When user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "Welch"
    When user moves the pointer away from box plot viewer
    And user selects "Student" in method choice input in box plot viewer
    Then "Method" property of box plot viewer should be "Student"
    When user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "Student"
    When user moves the pointer away from box plot viewer
    And user selects "Welch" in method choice input in box plot viewer
    Then "Method" property of box plot viewer should be "Welch"

  Scenario: Three groups, a control group and its comparisons table
    When user sets "Category 1" property of box plot viewer to "RACE"
    Then box plot viewer should have repainted
    And box plot viewer should have a "category Asian" area
    When user hovers over the "p value" area of box plot viewer
    Then tooltip should contain text "ANOVA"
    When user moves the pointer away from box plot viewer
    And user hovers over box plot viewer
    And user selects "Caucasian" in control-group choice input in box plot viewer
    Then "Control Comparisons" property of box plot viewer should be "true"
    And "Control Group" property of box plot viewer should be "Caucasian"
    And box plot viewer should have a "p value of Asian" area
    And box plot viewer should have a "control band" area
    When user picks "Add Control Comparisons Table" from the context menu of the "group comparison" area of box plot viewer
    Then table "Control Comparisons: AGE by RACE vs Caucasian" should be open
    And table "Control Comparisons: AGE by RACE vs Caucasian" should have 3 rows
    And table "Control Comparisons: AGE by RACE vs Caucasian" should have columns "Conclusion, Group, n, Mean, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g"
    And table "Control Comparisons: AGE by RACE vs Caucasian" should have no missing values in "p (adj)" column
    When user switches to the "demog-1000" table view
    And user clicks on the "p value of Asian" area of box plot viewer
    Then Results section in context panel should be present
    And Statistics section in context panel should contain text "Asian"

  Scenario: Two-way ANOVA
    When user sets properties of box plot viewer:
      | Control Comparisons | false |
      | Control Group       |       |
      | Category 2          | SEX   |
    And user hovers over box plot viewer
    Then baseline choice input in box plot viewer should be visible
    And box plot viewer should have a "RACE effect" area
    When user picks "Add Two-Way ANOVA Table" from the context menu of the "group comparison" area of box plot viewer
    Then table "Two-Way ANOVA: AGE by RACE, SEX" should be open
    And table "Two-Way ANOVA: AGE by RACE, SEX" should have 5 rows
    And table "Two-Way ANOVA: AGE by RACE, SEX" should have columns "Conclusion, Source of variance, SS, DF, MS, F, p-value"
    And table "Two-Way ANOVA: AGE by RACE, SEX" should have no missing values in "SS" column
    When user switches to the "demog-1000" table view

  Scenario: Closing the comparison
    When user hovers over the "group comparison" area of box plot viewer
    And user clicks on close-group-stats icon in box plot viewer
    Then "Show Group Comparison" property of box plot viewer should be "false"
    And "Show P Value" property of box plot viewer should be "true"
    And method choice input in box plot viewer should be hidden
    And box plot viewer should not have a "group comparison" area
    And box plot viewer should have a "p value" area

  Scenario: A covariate adjusts the value axis
    Given user sets properties of box plot viewer:
      | Category 2            |        |
      | Control Comparisons   | false  |
      | Control Group         |        |
      | Category 1            | SEX    |
      | Value                 | WEIGHT |
      | Show Group Comparison | true   |
    When user sets "Adjust By" property of box plot viewer to "HEIGHT"
    Then "Adjustment" property of box plot viewer should be "regressOut"
    When user hovers over box plot viewer
    Then "Adjust by" column input in box plot viewer should be visible
    And "Adjust by" column input in box plot viewer should contain text "HEIGHT"
    And Value column input in box plot viewer should have text "WEIGHT"
    And adjustment choice input in box plot viewer should have the value "Regress-out"
    When user selects "Ratio" in adjustment choice input in box plot viewer
    Then "Adjustment" property of box plot viewer should be "ratio"

  Scenario: ANCOVA against a control group
    When user sets "Category 1" property of box plot viewer to "RACE"
    And user hovers over box plot viewer
    And user selects "Caucasian" in control-group choice input in box plot viewer
    And user selects "ANCOVA" in method choice input in box plot viewer
    Then "Method" property of box plot viewer should be "ANCOVA"
    And "Adjustment" property of box plot viewer should be "ratio"
    And adjustment choice input in box plot viewer should be hidden
    When user hovers over box plot viewer
    Then "Adjust by" column input in box plot viewer should be visible
    And "Adjust by" column input in box plot viewer should contain text "HEIGHT"
    When user picks "Add ANCOVA Table" from the context menu of the "group comparison" area of box plot viewer
    Then table "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian" should be open
    And table "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian" should have 4 rows
    And table "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian" should have columns "Conclusion, Group, n, Raw mean, Adjusted mean, SE, Adj. diff, p-value, Hedges' g"
    And table "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian" should have no missing values in "Adjusted mean" column
    When user switches to the "demog-1000" table view

  Scenario: The matched baseline and the Simpson's paradox cue
    When user sets "Category 2" property of box plot viewer to "SEX"
    And user hovers over box plot viewer
    And user selects "Matched · per stratum" in baseline choice input in box plot viewer
    Then "Baseline Mode" property of box plot viewer should be "matched"
    And box plot viewer should have a "control band F" area
    And box plot viewer should have a "control band M" area
    And box plot viewer should not have a "control band" area
    When user adds a calculated column "SIMPSON_STRAT" with formula "if(Mod(Round(${HEIGHT} * 1000), 2) == 0, \"A\", \"B\")"
    And user adds a calculated column "SIMPSON_VAL" with formula "if(${SEX} == \"M\", 0, if(Mod(Round(${HEIGHT} * 1000), 2) == 0, 2.5, -2.5)) + (Mod(Round(${WEIGHT} * 137), 600) / 30)"
    And user sets properties of box plot viewer:
      | Category 1 | SEX           |
      | Category 2 | SIMPSON_STRAT |
      | Value      | SIMPSON_VAL   |
    And user hovers over box plot viewer
    And user selects "M" in control-group choice input in box plot viewer
    Then box plot viewer should have a "control band A" area
    And simpson-warning icon in box plot viewer should be visible
    When user hovers over simpson-warning icon in box plot viewer
    Then tooltip should contain text "Pooling cancels opposite within-stratum trends"
    When user sets properties of box plot viewer:
      | Value      | WEIGHT |
      | Category 1 | RACE   |
      | Category 2 | SEX    |
    And user removes "SIMPSON_STRAT" column
    And user removes "SIMPSON_VAL" column
    And user sets "Adjust By" property of box plot viewer to ""
    Then "Adjust By" property of box plot viewer should be ""
    And "Adjustment" property of box plot viewer should be ""
    And "Method" property of box plot viewer should be ""
