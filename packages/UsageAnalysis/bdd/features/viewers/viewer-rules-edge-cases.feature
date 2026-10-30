@viewers @edge-cases
Feature: Viewer edge cases shared by every viewer: no rows, a column removed or renamed, a layout with a missing column, a filter left by a closed viewer, values out of bounds, empty categories
  One rule per scenario outline, one Examples row per viewer it applies to. A row that breaks its
  rule is a known failure in an Examples block of its own, tagged with its ticket. demog-1000: 1000
  rows, none with AGE between 500 and 600. The Charts viewers' rows are in packages/Charts/bdd.
  - No rows: with every row filtered out a viewer shows no row, draws nothing, reports no error and
    shows no message; the rows come back with the filter.
  - A column removed: a column a viewer uses for an encoding (color, size, lanes, category) leaves
    the viewer's settings and nothing is logged. A viewer may pick another column in its place (the
    box plot picks another category). The bar chart's Stack is in bar-chart-edge-cases.feature.
  - A column renamed: every setting that named the column follows the new name, a Filter formula
    included.
  - A layout naming a missing column: the rest of the layout is applied, the viewer draws without
    the column and nothing is logged. The setting may keep the missing column's name.
  - A filter left by a viewer: the rows a viewer's click filtered out come back when the viewer is
    closed, or when the column of the clicked category is removed.
  - Out of bounds: a value typed in the context panel outside the property's declared bounds is not
    taken; the property keeps its value. The panel's number input marks such a value invalid
    (`number_input.dart`, `addMinMaxValidators`), but the input writes every change into the
    property (`InputBase.bind`).
  - Empty categories: an empty value is a category named `(empty)`, as in the tree.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset

  Scenario Outline: With every row filtered out a <viewer> shows no row and no message, and the rows come back
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user filters rows where "AGE" is between 500 and 600
    Then 0 rows should pass the filter
    And the "rows shown" reading of <viewer> viewer should be 0
    And the "<drawn>" reading of <viewer> viewer should be 0
    And <viewer> viewer should report no error
    When user resets the filter
    Then the "rows shown" reading of <viewer> viewer should be 1000
    And the "<drawn>" reading of <viewer> viewer should be at least 1
    And no errors should have been logged

    Examples:
      | viewer          | property             | value               | drawn           |
      | bar chart       | splitColumnName      | RACE                | rows shown      |
      | box plot        | valueColumnName      | AGE                 | rows shown      |
      | heat map        | rowSource            | Filtered            | rows shown      |
      | histogram       | valueColumnName      | AGE                 | rows shown      |
      | matrix plot     | rowSource            | Filtered            | rows shown      |
      | network diagram | node1ColumnName      | RACE                | rows shown      |
      | pc plot         | columnNames          | AGE, HEIGHT, WEIGHT | lines drawn     |
      | pie chart       | categoryColumnName   | RACE                | slices          |
      | pivot table     | rowSource            | Filtered            | aggregated rows |
      | scatter plot    | xColumnName          | AGE                 | rows shown      |
      | tile viewer     | lanesColumnName      | RACE                | tiles           |
      | tree map        | splitByColumnNames   | RACE, SEX           | leaves          |
      | trellis plot    | yColumnNames         | RACE                | rows shown      |

  Scenario Outline: The <property> column removed from the table leaves the settings of the <viewer>
    When user adds a <viewer> viewer with:
      | <property> | <column> |
    And user removes "<column>" column
    Then no errors should have been logged
    And "<property>" property of <viewer> viewer should not be "<column>"
    And <viewer> viewer should report no error

    Examples:
      | viewer       | property            | column |
      | scatter plot | colorColumnName     | RACE   |
      | scatter plot | sizeColumnName      | HEIGHT |
      | tile viewer  | lanesColumnName     | RACE   |
      | tree map     | colorColumnName     | AGE    |
      | tree map     | sizeColumnName      | WEIGHT |
      | box plot     | category1ColumnName | SEX    |

  Scenario Outline: The <column> column renamed is followed by the <property> setting of the <viewer>
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user renames "<column>" column to "<new name>"
    Then "<property>" property of <viewer> viewer should be "<renamed>"
    And <viewer> viewer should report no error
    And no errors should have been logged

    Examples:
      | viewer       | property             | value               | column | new name  | renamed              |
      | scatter plot | xColumnName          | AGE                 | AGE    | AGE_YEARS | AGE_YEARS            |
      | scatter plot | colorColumnName      | RACE                | RACE   | ORIGIN    | ORIGIN               |
      | pie chart    | categoryColumnName   | RACE                | RACE   | ORIGIN    | ORIGIN               |
      | pc plot      | columnNames          | AGE, HEIGHT, WEIGHT | HEIGHT | TALL      | AGE, TALL, WEIGHT    |
      | tile viewer  | lanesColumnName      | RACE                | RACE   | ORIGIN    | ORIGIN               |
      | tree map     | splitByColumnNames   | RACE, SEX           | SEX    | GENDER    | RACE, GENDER         |
      | trellis plot | yColumnNames         | RACE                | RACE   | ORIGIN    | ORIGIN               |

    @known-failure:78 @realizes:core#1213
    Examples: A Filter formula keeps the old name: the rename handler renames column settings only (`DataFrameViewer._renameColumnInLook`)
      | viewer       | property | value       | column | new name | renamed       |
      | scatter plot | filter   | ${AGE} > 44 | AGE    | YEARS    | ${YEARS} > 44 |

  Scenario Outline: A layout naming a <column> column the table no longer has is applied to the <viewer>, which draws without it
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user saves the layout of the current table view
    And user removes "<column>" column
    And user loads the saved layout
    Then <viewer> viewer should be visible
    And "<property>" property of <viewer> viewer should be "<left>"
    And <viewer> viewer should report no error
    And no errors should have been logged

    Examples: The setting keeps the missing column's name
      | viewer       | property             | value               | column | left                |
      | scatter plot | colorColumnName      | RACE                | RACE   | RACE                |
      | pc plot      | columnNames          | AGE, HEIGHT, WEIGHT | HEIGHT | AGE, HEIGHT, WEIGHT |
      | tree map     | splitByColumnNames   | RACE, SEX           | SEX    | RACE, SEX           |

    @known-failure:106 @realizes:public#4452
    Examples: The forms viewer takes the removed column's name back into its fields and throws a TypeError on every render
      | viewer       | property             | value               | column | left                |
      | forms        | fields               | AGE, HEIGHT, SEX    | HEIGHT | AGE, HEIGHT, SEX    |

  Scenario Outline: Closing a <viewer> gives back the rows its click filtered out
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user sets "onClick" property of <viewer> viewer to "Filter"
    And user clicks on the "<area>" area of <viewer> viewer
    Then fewer than 1000 rows should pass the filter
    When user clicks on close icon of <viewer> viewer
    Then the open tableview should have 0 <viewer> viewers
    And all rows should pass the filter
    And no errors should have been logged

    Examples:
      | viewer    | property        | value | area      |
      | bar chart | splitColumnName | RACE  | bar Asian |

    @known-failure:127 @realizes:core#233
    Examples: The calendar leaves the table filtered until something else asks for a refilter
      | viewer   | property | value  | area           |
      | calendar | onClick  | Filter | month 1990-01  |

  Scenario Outline: Removing the column of the category a <viewer> click filtered by gives the rows back
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user sets "onClick" property of <viewer> viewer to "Filter"
    And user clicks on the "<area>" area of <viewer> viewer
    Then fewer than 1000 rows should pass the filter
    When user removes "<column>" column
    Then all rows should pass the filter
    And no errors should have been logged

    @known-failure:146 @realizes:core#1174
    Examples: The removal does not refilter the table; the rows come back on the next refilter
      | viewer    | property             | value              | area                            | column |
      | bar chart | splitColumnName      | RACE               | bar Asian                       | RACE   |

  Scenario Outline: <caption> typed in the context panel outside its declared bounds is not taken by the <viewer>
    When user adds a <viewer> viewer
    And user clicks on settings icon of <viewer> viewer
    Then context panel should be visible
    Given "<category>" category in context panel is expanded
    When user enters "<typed>" into "<caption>" property in context panel
    Then "<caption>" property of <viewer> viewer should be "<kept>"
    And no errors should have been logged

    @known-failure:160 @realizes:core#1200
    Examples: The value is written into the property although the input marks it invalid
      | viewer       | category | caption            | typed | kept |
      | histogram    | X Axis   | Bins               | 201   | 20   |
      | pie chart    | Misc     | Max Radius         | 600   | 150  |
      | trellis plot | Misc     | Tiles Per Row      | 11    | 4    |

  Scenario Outline: <property> set through the API outside its declared bounds is not taken by the <viewer>
    When user adds a <viewer> viewer
    And user sets "<property>" property of <viewer> viewer to "<set>"
    Then "<property>" property of <viewer> viewer should be "<kept>"
    And no errors should have been logged

    @known-failure:173 @realizes:core#1200
    Examples: The value is written into the property
      | viewer       | property         | set | kept |
      | histogram    | bins             | 201 | 20   |
      | histogram    | bins             | 0   | 20   |
      | pie chart    | maxRadius        | 600 | 150  |
      | trellis plot | tilesPerRow      | 0   | 4    |

  Scenario Outline: An empty value is the (empty) category of the <viewer>
    Given user opens a table "with-empty" with:
      | site | arm |
      | A    | x   |
      | A    |     |
      | B    | y   |
    When user adds a <viewer> viewer with:
      | <property> | arm |
    Then no errors should have been logged
    And <viewer> viewer should report no error
    And the "<shown>" reading of <viewer> viewer should be 3
    And the "<reading>" reading of <viewer> viewer should contain "(empty)"

    @known-failure:195 @realizes:core#1217
    Examples: The tile viewer shows the empty value as a lane with no name ("x, y, ")
      | viewer      | property        | shown | reading    |
      | tile viewer | lanesColumnName | tiles | lane names |
