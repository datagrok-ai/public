@journey @viewers @realizes:viewers.word-cloud
Feature: Word cloud counts, the category gate and the viewer filter
  What the cloud is: one word per distinct value of a string column, sized by how many rows carry
  it. Every claim here is the viewer's own count — `rows of word "Caucasian"` — where the spec this
  replaces read `chart.getOption().series[0].data`, a private echarts structure, and then re-counted
  the column in JavaScript to check it against itself. A word is pointed at by name (`word "F"`),
  not by scanning the canvas in 50×30 pixel buckets and hovering the six densest until a tooltip
  appeared.
  Two facts of the fixture are asserted in the Background because everything else rests on them:
  demog-1000's first countable string column is SEX (F 553 / M 447 — USUBJID is skipped, it has
  1000 distinct values), and RACE is Caucasian 896 / Other 62 / Black 27 / Asian 15.
  The 500-category gate is a message, not an empty cloud: `render` returns before creating the
  chart, so the viewer reports `error` and NO word areas and NO cloud readings, while `this.chart`
  still holds the previous cloud's geometry. That trap is the point of the scenario that checks it.
  The filtering scenario is rewritten. The old one compared the word NAMES before and after and
  declared the cloud alive; the viewer counts over its own filter, so the names are exactly what a
  filter that leaves one row of each category cannot change — the counts are what must move.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a word cloud viewer
    Then word cloud viewer should be visible
    And 1000 rows should pass the filter
    And the "column" reading of word cloud viewer should be "SEX"
    And word cloud viewer should be painted

  Scenario: The cloud draws one word per distinct value with that value's row count
    Then the "words" reading of word cloud viewer should be 2
    And the "word names" reading of word cloud viewer should be "F, M"
    And the "rows of word \"F\"" reading of word cloud viewer should be 553
    And the "rows of word \"M\"" reading of word cloud viewer should be 447
    And the "rows shown" reading of word cloud viewer should be 1000
    And word cloud viewer should have a "word \"F\"" area
    And word cloud viewer should have a "word \"M\"" area
    And word cloud viewer should report no error
    And no errors should have been logged

  Scenario: Column RACE draws one word per race with its row count
    When user sets "wordColumnName" property of word cloud viewer to "RACE"
    Then the "words" reading of word cloud viewer should be 4
    And the "word names" reading of word cloud viewer should contain "Caucasian"
    And the "word names" reading of word cloud viewer should contain "Asian"
    And the "rows of word \"Caucasian\"" reading of word cloud viewer should be 896
    And the "rows of word \"Other\"" reading of word cloud viewer should be 62
    And the "rows of word \"Black\"" reading of word cloud viewer should be 27
    And the "rows of word \"Asian\"" reading of word cloud viewer should be 15
    And word cloud viewer should have repainted
    And no errors should have been logged

  Scenario: The word of the bigger count gets the bigger box
    Rotation is pinned to zero so that a box's height is the word's text size and nothing else: at
    the default -30..30 each word is drawn at a random angle, and the axis-aligned box a rotated
    word reports mixes its width into its height.
    Given user sets properties of word cloud viewer:
      | minRotationDegree | 0 |
      | maxRotationDegree | 0 |
    Then the "word \"Caucasian\"" area of word cloud viewer should be taller than the "word \"Other\"" area
    And the "word \"Other\"" area of word cloud viewer should be taller than the "word \"Asian\"" area
    And the "word \"Caucasian\"" area of word cloud viewer should be at least 60 pixels tall
    When user sets properties of word cloud viewer:
      | minRotationDegree | -30 |
      | maxRotationDegree | 30  |
    Then the "words" reading of word cloud viewer should be 4
    And no errors should have been logged

  Scenario: A column with more than 500 categories reports the gate and shows no cloud at all
    When user sets "wordColumnName" property of word cloud viewer to "USUBJID"
    Then word cloud viewer should report the error "The Word cloud viewer requires categorical column with 500 or fewer unique categories"
    And word cloud viewer should not have a "word \"Caucasian\"" area
    And word cloud viewer should not have a "view" area
    And word cloud viewer should not report a "words" reading
    And word cloud viewer should not report a "word names" reading
    And word cloud viewer should not report a "rows of word \"Caucasian\"" reading
    And the "column" reading of word cloud viewer should be "USUBJID"
    And the "rows shown" reading of word cloud viewer should be 1000
    When user sets "wordColumnName" property of word cloud viewer to "RACE"
    Then word cloud viewer should report no error
    And the "words" reading of word cloud viewer should be 4
    And word cloud viewer should have a "word \"Caucasian\"" area
    And no errors should have been logged

  Scenario: Hovering a word shows how many rows carry it
    The rotation range is pinned to zero for this one: a word is drawn at a random angle inside the
    range, and the area a rotated word reports is the axis-aligned box around it, whose centre can
    fall beside the glyphs for a short word. Upright, the box is the word.
    Given user sets properties of word cloud viewer:
      | minRotationDegree | 0 |
      | maxRotationDegree | 0 |
    When user hovers over the "word \"Caucasian\"" area of word cloud viewer
    Then exactly one tooltip should be shown
    And tooltip should contain text "896 rows"
    When user moves the pointer away from word cloud viewer
    And user hovers over the "word \"Other\"" area of word cloud viewer
    Then exactly one tooltip should be shown
    And tooltip should contain text "62 rows"
    When user moves the pointer away from word cloud viewer
    And user sets properties of word cloud viewer:
      | minRotationDegree | -30 |
      | maxRotationDegree | 30  |
    Then the "words" reading of word cloud viewer should be 4
    And no errors should have been logged

  Scenario: Clicking a word selects exactly that word's rows
    Given user clears the row selection
    When user clicks on the "word \"Other\"" area of word cloud viewer
    Then 62 rows should be selected
    And only rows where "RACE" is "Other" should be selected
    When user clicks on the "word \"Asian\"" area of word cloud viewer
    Then 15 rows should be selected
    And only rows where "RACE" is "Asian" should be selected
    When user clears the row selection
    Then no rows should be selected
    And no errors should have been logged

  Scenario: A table filter moves the counts and leaves the names alone
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of word cloud viewer should be 447
    And the "rows of word \"Caucasian\"" reading of word cloud viewer should be 416
    And the "rows of word \"Other\"" reading of word cloud viewer should be 14
    And the "words" reading of word cloud viewer should be 4
    And the "word names" reading of word cloud viewer should contain "Caucasian"
    And word cloud viewer should report no error
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows of word \"Caucasian\"" reading of word cloud viewer should be 896
    And the "rows shown" reading of word cloud viewer should be 1000
    And no errors should have been logged

  Scenario: The title bar closes the cloud
    When user clicks on close icon of word cloud viewer
    Then word cloud viewer should be absent
    And the open tableview should have 0 word cloud viewers
    And no errors should have been logged
