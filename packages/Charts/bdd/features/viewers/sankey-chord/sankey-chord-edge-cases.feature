@viewers @edge-cases @realizes:charts.viewer.sankey @realizes:charts.viewer.chord
Feature: Sankey and Chord edge cases — Source equal to Target, too many categories, a removed column
  The Sankey and the Chord on settings and data changes their other features do not reach. A
  Sankey whose Source and Target are the same column links every category to itself, and the
  viewer says "The graph contains cycles. Please remove circular dependencies." (`sankey.ts`, the
  cycle check before drawing). A Chord whose From and To columns hold more than 1500 category
  pairs says "Too many categories to render." (`chord-viewer.ts`, `maxLinkNumber`): USUBJID holds
  1000 categories on demog-1000, and with SEX that is 2000 pairs. A Source column removed from the
  table while the Sankey shows it, and a From column removed while the Chord shows it, should leave
  the viewers without an error; both log "No such column" (public#4420, public#4421). demog-1000; the Sankey
  opens on SEX, RACE and AGE, the Chord on SEX and RACE.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset

  Scenario: A Sankey whose Source and Target are the same column says the graph has cycles
    When user adds a sankey viewer
    Then "Source" property of sankey viewer should be "SEX"
    When user sets "Target" property of sankey viewer to "SEX"
    Then "The graph contains cycles. Please remove circular dependencies." text in sankey viewer should be visible
    And no errors should have been logged
    When user sets "Target" property of sankey viewer to "RACE"
    Then "Caucasian" text in sankey viewer should be visible
    And "The graph contains cycles. Please remove circular dependencies." text in sankey viewer should be absent
    And no errors should have been logged

  Scenario: A Chord over more than 1500 category pairs says there are too many categories
    When user adds a chord viewer
    Then "From" property of chord viewer should be "SEX"
    When user sets "To" property of chord viewer to "USUBJID"
    Then "Too many categories to render." text in chord viewer should be visible
    And no errors should have been logged
    When user sets "To" property of chord viewer to "RACE"
    Then "Asian" text in chord viewer should be visible
    And "Too many categories to render." text in chord viewer should be absent
    And no errors should have been logged

  @realizes:public#4420
  Scenario: The Sankey's Source column removed from the table raises no error
    Guards public#4420 (fixed): the Sankey kept reading the removed column and logged "No such
    column: SEX" (`sankey.ts`, `getCol` in `prepareData`).
    When user adds a sankey viewer
    Then "Source" property of sankey viewer should be "SEX"
    When user removes "SEX" column
    Then the table should not have a column "SEX"
    When user filters rows where "AGE" is between 30 and 50
    Then 494 rows should pass the filter
    And no errors should have been logged

  @realizes:public#4421
  Scenario: The Chord's From column removed from the table raises no error
    Guards public#4421 (fixed): the Chord kept reading the removed column and logged "No such
    column: SEX" (`chord-viewer.ts`, `getCol` in `draw`).
    When user adds a chord viewer
    Then "From" property of chord viewer should be "SEX"
    When user removes "SEX" column
    Then the table should not have a column "SEX"
    When user filters rows where "AGE" is between 30 and 50
    Then 494 rows should pass the filter
    And no errors should have been logged

  @known-failure:73 @realizes:public#4462
  Scenario: A Sankey whose Value column holds negative values refuses to draw and says so, as the pie does
    The pie chart refuses a negative angle column with "sum(NEG) contains negative values"
    (`pie_chart_core.dart`, `_aggregate`); the Sankey should refuse the same way. The claim reads
    the words the two messages share.
    When user adds a calculated column "NEG" with formula "-${AGE}"
    And user adds a sankey viewer
    Then "Value" property of sankey viewer should be "AGE"
    When user sets "Value" property of sankey viewer to "NEG"
    Then sankey viewer should contain the text "contains negative values"
    And no errors should have been logged
    When user sets "Value" property of sankey viewer to "AGE"
    Then "Caucasian" text in sankey viewer should be visible
    And sankey viewer should not contain the text "contains negative values"
    And no errors should have been logged

  @known-failure:88 @realizes:public#4462
  Scenario: A Sankey whose Value column holds only zeros refuses to draw and says so, as the pie does
    The pie chart refuses an angle column of zeros with "sum(ZERO): all values are 0"
    (`pie_chart_core.dart`, `_aggregate`); the Sankey should refuse the same way.
    When user adds a calculated column "ZERO" with formula "0"
    And user adds a sankey viewer
    Then "Value" property of sankey viewer should be "AGE"
    When user sets "Value" property of sankey viewer to "ZERO"
    Then sankey viewer should contain the text "all values are 0"
    And no errors should have been logged
    When user sets "Value" property of sankey viewer to "AGE"
    Then "Caucasian" text in sankey viewer should be visible
    And sankey viewer should not contain the text "all values are 0"
    And no errors should have been logged
