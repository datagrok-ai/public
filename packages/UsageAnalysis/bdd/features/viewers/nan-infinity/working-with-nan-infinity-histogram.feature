@journey @viewers @realizes:viewers.histogram
Feature: A histogram over a column that holds a NaN and an Infinity
  `working-with-nan-infinity.md` asks for viewers over the HEIGHT and WEIGHT columns once a NaN and
  an Infinity are written into them; `working-with-nan-infinity.feature` walks the scatter plot, and
  this is the histogram. An Infinity in WEIGHT must not stretch the value axis to a number it cannot
  draw: the axis keeps its finite ends, the bins stay where they were, and the viewer reports no
  error; a NaN in HEIGHT leaves that axis finite as well. One journey on demog-1000, where every
  WEIGHT is set and 128 HEIGHTs are blank; the last scenario puts both values back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | valueColumnName | WEIGHT |
    Then histogram viewer should report no error
    And the "bins shown" reading of histogram viewer should be 20

  Scenario: An Infinity in WEIGHT keeps the value axis finite and the bins as they were
    When user remembers the "axis max" reading of histogram viewer
    And user takes a snapshot of histogram viewer
    And user sets "WEIGHT" column in row 2 to "Infinity"
    Then histogram viewer should have repainted
    And the "axis max" reading of histogram viewer should be a finite number
    And the "axis max" reading of histogram viewer should be as remembered
    And the "bins shown" reading of histogram viewer should be 20
    And histogram viewer should report no error
    And histogram viewer should be painted
    And no errors should have been logged

  Scenario: A NaN in HEIGHT leaves the axis finite
    When user sets "valueColumnName" property of histogram viewer to "HEIGHT"
    Then the "axis max" reading of histogram viewer should differ from before
    When user remembers the "axis max" reading of histogram viewer
    And user takes a snapshot of histogram viewer
    And user sets "HEIGHT" column in row 1 to "NaN"
    Then histogram viewer should have repainted
    And the "axis min" reading of histogram viewer should be a finite number
    And the "axis max" reading of histogram viewer should be as remembered
    And histogram viewer should report no error
    And no errors should have been logged

  Scenario: The two values put back leave the axis as it was
    When user sets "HEIGHT" column in row 1 to "174.705"
    Then the "axis max" reading of histogram viewer should be as remembered
    When user sets "WEIGHT" column in row 2 to "64"
    And user sets "valueColumnName" property of histogram viewer to "WEIGHT"
    Then the "axis max" reading of histogram viewer should be 165
    And histogram viewer should report no error
    And no errors should have been logged
