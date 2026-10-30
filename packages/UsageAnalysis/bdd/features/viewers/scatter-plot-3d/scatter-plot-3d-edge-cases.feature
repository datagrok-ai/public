@viewers @edge-cases @realizes:viewers.scatter-plot3d
Feature: 3D scatter plot edge cases: no rows, columns removed or renamed, a table swap
  Guards for the changes under an open 3D scatter plot, claimed on what the viewer reports (the
  WebGL picture has no readable pixels): `rows shown`, the `scene signature` of the frame, the
  column settings and the on-viewer column selectors. A filter that removes every row leaves the
  plot with nothing to draw and no error, and the frame comes back with the rows. A renamed column
  is followed by every column setting that named it (the viewer base renames column settings on the
  table's rename event). Three defects of core#1027 are pinned as known failures: removing the Size
  column throws in the scene update, removing the Color column resets the X, Y and Z the user chose,
  and a table swap and back takes the four column selectors off the plot. Each tag goes when
  core#1027 is fixed. Data: demog-1000, X = AGE, Y = WEIGHT, Z = HEIGHT (872 rows with all three).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a 3d scatter plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
      | zColumnName | HEIGHT |
    Then the "rows shown" reading of 3d scatter plot viewer should be 872

  Scenario: A filter that removes every row leaves the plot empty without an error
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of 3d scatter plot viewer should be 0
    And the "scene signature" reading of 3d scatter plot viewer should not be as remembered
    When user resets the filter
    Then the "rows shown" reading of 3d scatter plot viewer should be 872
    And the "scene signature" reading of 3d scatter plot viewer should be as remembered
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Renaming the axis columns is followed by the plot's settings
    When user renames "AGE" column to "AGE_YEARS"
    And user renames "HEIGHT" column to "HEIGHT_CM"
    Then "xColumnName" property of 3d scatter plot viewer should be "AGE_YEARS"
    And "zColumnName" property of 3d scatter plot viewer should be "HEIGHT_CM"
    And the "rows shown" reading of 3d scatter plot viewer should be 872
    And no errors should have been logged

  @realizes:core#1027
  Scenario: Removing the Size column keeps the plot working without an error
    Guards core#1027 (fixed): removing the size column threw in updateAll, and the viewer then
    ignored every later property change.
    When user sets "Size" property of 3d scatter plot viewer to "WEIGHT"
    And user removes "WEIGHT" column
    Then no errors should have been logged
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user sets "Marker Type" property of 3d scatter plot viewer to "box"
    Then the "scene signature" reading of 3d scatter plot viewer should not be as remembered

  @realizes:core#1027
  Scenario: Removing the Color column keeps the chosen X, Y and Z
    Guards core#1027 (fixed): removing only the color column reset the X, Y and Z the user chose.

    When user sets "Color" property of 3d scatter plot viewer to "RACE"
    And user removes "RACE" column
    Then "xColumnName" property of 3d scatter plot viewer should be "AGE"
    And "yColumnName" property of 3d scatter plot viewer should be "WEIGHT"
    And "zColumnName" property of 3d scatter plot viewer should be "HEIGHT"
    And no errors should have been logged

  @realizes:core#1027
  Scenario: A table swap and back keeps the on-viewer column selectors
    Guards core#1027 (fixed): a table swap removed the four column selectors from the plot.

    Given user opens cars dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of 3d scatter plot viewer to "cars"
    Then 3d scatter plot viewer should be bound to table "cars"
    When user sets "table" property of 3d scatter plot viewer to "demog-1000"
    Then 3d scatter plot viewer should be bound to table "demog-1000"
    And X column input in 3d scatter plot viewer should be visible
    And Y column input in 3d scatter plot viewer should be visible
    And Z column input in 3d scatter plot viewer should be visible
    And Color column input in 3d scatter plot viewer should be visible
