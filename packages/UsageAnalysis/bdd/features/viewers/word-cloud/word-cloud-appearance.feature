@journey @viewers @realizes:viewers.word-cloud
Feature: Word cloud text size, rotation, shape and font
  The four appearance properties, claimed on the boxes the viewer reports rather than on the canvas
  as a whole. Three of the scenarios this replaces asserted only that "more than 100 pixels
  changed" after the property was set — true of any two different layouts of the same words, and so
  true whether or not the property did anything. Equal min and max text size is claimed by the word
  boxes becoming one height; a pinned rotation is claimed by the boxes turning; and Shape is
  claimed by the cloud keeping every word and its count, because the mask decides what still fits.
  Rotation is pinned before every box claim. The layout draws each word at a random angle inside
  the rotation range, so two renders of the SAME cloud have different boxes — a box read under the
  default -30..30 says nothing repeatable, which is the reason the scenarios below set the range to
  a single value first and restore it after.
  The Background pins RACE: Caucasian 896, Other 62, Black 27, Asian 15 on demog-1000 — four words
  whose sizes are far enough apart for a size claim to mean something.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a word cloud viewer with:
      | wordColumnName    | RACE |
      | minRotationDegree | 0    |
      | maxRotationDegree | 0    |
    Then word cloud viewer should be visible
    And the "words" reading of word cloud viewer should be 4
    And the "rows of word \"Caucasian\"" reading of word cloud viewer should be 896
    And word cloud viewer should be painted

  Scenario: Text size spreads the boxes by count, and equal min and max flattens them
    Then the "word \"Caucasian\"" area of word cloud viewer should be taller than the "word \"Other\"" area
    And the "word \"Other\"" area of word cloud viewer should be taller than the "word \"Black\"" area
    When user sets properties of word cloud viewer:
      | minTextSize | 20 |
      | maxTextSize | 20 |
    Then the "word \"Caucasian\"" and "word \"Black\"" areas of word cloud viewer should be the same height
    And the "word \"Caucasian\"" and "word \"Other\"" areas of word cloud viewer should be the same height
    And the "word \"Caucasian\"" and "word \"Asian\"" areas of word cloud viewer should be the same height
    And the "word \"Caucasian\"" area of word cloud viewer should be shorter than before
    And the "word \"Caucasian\"" area of word cloud viewer should be wider than the "word \"Black\"" area
    When user sets properties of word cloud viewer:
      | minTextSize | 14  |
      | maxTextSize | 100 |
    Then the "word \"Caucasian\"" area of word cloud viewer should be taller than before
    And the "word \"Caucasian\"" area of word cloud viewer should be taller than the "word \"Black\"" area
    And no errors should have been logged

  Scenario: A rotation of 90 degrees turns every box on its side
    When user sets properties of word cloud viewer:
      | minTextSize | 20 |
      | maxTextSize | 20 |
    Then the "word \"Caucasian\"" area of word cloud viewer should be wider than the "word \"Black\"" area
    When user sets properties of word cloud viewer:
      | minRotationDegree | 90 |
      | maxRotationDegree | 90 |
    Then the "word \"Caucasian\"" area of word cloud viewer should be narrower than before
    And the "word \"Caucasian\"" area of word cloud viewer should be taller than before
    And the "word \"Caucasian\"" and "word \"Black\"" areas of word cloud viewer should be the same width
    And the "word \"Caucasian\"" area of word cloud viewer should be taller than the "word \"Black\"" area
    When user sets properties of word cloud viewer:
      | minRotationDegree | 0   |
      | maxRotationDegree | 0   |
      | minTextSize       | 14  |
      | maxTextSize       | 100 |
    Then the "word \"Caucasian\"" area of word cloud viewer should be wider than the "word \"Black\"" area
    And no errors should have been logged

  Scenario: Every shape keeps all four words and their counts
    Then "shape" property of word cloud viewer should be "circle"
    When user sets "shape" property of word cloud viewer to "diamond"
    Then the "words" reading of word cloud viewer should be 4
    And the "rows of word \"Caucasian\"" reading of word cloud viewer should be 896
    And word cloud viewer should have a "word \"Asian\"" area
    And word cloud viewer should have repainted
    When user sets "shape" property of word cloud viewer to "star"
    Then the "words" reading of word cloud viewer should be 4
    And the "rows of word \"Asian\"" reading of word cloud viewer should be 15
    And word cloud viewer should have a "word \"Black\"" area
    When user sets "shape" property of word cloud viewer to "pentagon"
    Then the "words" reading of word cloud viewer should be 4
    And word cloud viewer should have a "word \"Other\"" area
    When user sets "shape" property of word cloud viewer to "triangle"
    Then the "words" reading of word cloud viewer should be 4
    And word cloud viewer should have a "word \"Caucasian\"" area
    When user sets "shape" property of word cloud viewer to "circle"
    Then the "words" reading of word cloud viewer should be 4
    And word cloud viewer should report no error
    And no errors should have been logged

  Scenario: The font reading is the font the words were drawn with, not the property
    Then the "font" reading of word cloud viewer should be "bold sans-serif"
    When user sets "bold" property of word cloud viewer to "false"
    Then the "font" reading of word cloud viewer should be "normal sans-serif"
    When user sets "fontFamily" property of word cloud viewer to "monospace"
    Then the "font" reading of word cloud viewer should be "normal monospace"
    And the "words" reading of word cloud viewer should be 4
    When user sets properties of word cloud viewer:
      | bold       | true       |
      | fontFamily | sans-serif |
    Then the "font" reading of word cloud viewer should be "bold sans-serif"
    And no errors should have been logged
