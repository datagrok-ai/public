@viewers @edge-cases
Feature: Charts viewer edge cases shared by every viewer: no rows, a column removed or renamed, a layout with a missing column, a filter left by a closed viewer, values out of bounds, empty categories, a category limit
  The rules of `packages/UsageAnalysis/bdd/features/viewers/viewer-rules-edge-cases.feature`, with
  the rows of the Charts viewers (sunburst, tree, word cloud). One rule per scenario outline, one
  Examples row per viewer it applies to. A row that breaks its rule is a known failure in an
  Examples block of its own, tagged with its ticket. demog-1000: 1000 rows, none with AGE between
  500 and 600.
  - No rows: with every row filtered out a viewer shows no row, reports no error and shows no
    message; the rows come back with the filter.
  - A column removed: a column a viewer uses for an encoding (color, size) leaves the viewer's
    settings and nothing is logged. A viewer may pick another column in its place.
  - A column renamed: every setting that named the column follows the new name.
  - A layout naming a missing column: the rest of the layout is applied, the viewer draws without
    the column and nothing is logged. The setting may keep the missing column's name.
  - A filter left by a viewer: the rows a viewer's click filtered out come back when the viewer is
    closed, or when the column of the clicked category is removed.
  - Out of bounds: a value typed in the context panel outside the property's declared bounds is not
    taken; the property keeps its value.
  - Empty categories: an empty value is a category named `(empty)`, as in the tree.
  - A category limit: a column with more than 500 categories is refused whatever the filter. Today
    the tree and the sunburst check the limit only while at least 500 rows pass the filter
    (`tree-viewer.ts` and `sunburst-viewer.ts`, `filter.trueCount >= CATEGORIES_NUMBER`).

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
    And <viewer> viewer should report no error
    When user resets the filter
    Then the "rows shown" reading of <viewer> viewer should be 1000
    And no errors should have been logged

    Examples:
      | viewer     | property             | value              |
      | sunburst   | hierarchyColumnNames | SEX, RACE          |
      | tree       | hierarchyColumnNames | CONTROL, SEX, RACE |
      | word cloud | wordColumnName       | RACE               |

  Scenario Outline: The <property> column removed from the table leaves the settings of the <viewer>
    When user adds a <viewer> viewer with:
      | <property> | <column> |
    And user removes "<column>" column
    Then no errors should have been logged
    And "<property>" property of <viewer> viewer should not be "<column>"
    And <viewer> viewer should report no error

    @realizes:public#4423
    Examples: The tree drops the removed Size or Color column (public#4423, fixed)
      | viewer | property        | column |
      | tree   | sizeColumnName  | HEIGHT |
      | tree   | colorColumnName | AGE    |

  Scenario Outline: The <column> column renamed is followed by the <property> setting of the <viewer>
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user renames "<column>" column to "<new name>"
    Then "<property>" property of <viewer> viewer should be "<renamed>"
    And <viewer> viewer should report no error
    And no errors should have been logged

    Examples:
      | viewer     | property             | value              | column | new name | renamed              |
      | tree       | hierarchyColumnNames | CONTROL, SEX, RACE | RACE   | ORIGIN   | CONTROL, SEX, ORIGIN |
      | word cloud | wordColumnName       | RACE               | RACE   | ORIGIN   | ORIGIN               |
      | sunburst   | hierarchyColumnNames | SEX, RACE          | RACE   | ORIGIN   | SEX, ORIGIN          |

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
      | viewer   | property             | value              | column | left               |
      | sunburst | hierarchyColumnNames | SEX, RACE          | SEX    | SEX, RACE          |
      | tree     | hierarchyColumnNames | CONTROL, SEX, RACE | RACE   | CONTROL, SEX, RACE |

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

    @realizes:public#4425
    Examples: Closing the sunburst or the tree in Filter mode gives the rows back (public#4425, fixed)
      | viewer   | property             | value              | area                            |
      | sunburst | hierarchyColumnNames | SEX, RACE          | segment F                       |
      | tree     | hierarchyColumnNames | CONTROL, SEX, RACE | branch All \| false \| F \| Asian |

  Scenario Outline: Removing the column of the category a <viewer> click filtered by gives the rows back
    When user adds a <viewer> viewer with:
      | <property> | <value> |
    And user sets "onClick" property of <viewer> viewer to "Filter"
    And user clicks on the "<area>" area of <viewer> viewer
    Then fewer than 1000 rows should pass the filter
    When user removes "<column>" column
    Then all rows should pass the filter
    And no errors should have been logged

    @known-failure:114 @realizes:core#1174
    Examples: The removal does not refilter the table; the rows come back on the next refilter
      | viewer   | property             | value              | area                            | column |
      | tree     | hierarchyColumnNames | CONTROL, SEX, RACE | branch All \| false \| F \| Asian | RACE   |
      | sunburst | hierarchyColumnNames | SEX, RACE          | segment F                       | SEX    |

  Scenario Outline: <caption> typed in the context panel outside its declared bounds is not taken by the <viewer>
    When user adds a <viewer> viewer
    And user clicks on settings icon of <viewer> viewer
    Then context panel should be visible
    Given "<category>" category in context panel is expanded
    When user enters "<typed>" into "<caption>" property in context panel
    Then "<caption>" property of <viewer> viewer should be "<kept>"
    And no errors should have been logged

    @known-failure:129 @realizes:core#1200
    Examples: The value is written into the property although the input marks it invalid
      | viewer | category | caption            | typed | kept |
      | tree   | Style    | Initial Tree Depth | 6     | 3    |
      | tree   | Style    | Font Size          | 31    | 12   |

    @known-failure:129 @realizes:public#4424
    Examples: Rotation Step 0 is taken and the cloud goes blank
      | viewer     | category | caption       | typed | kept |
      | word cloud | Misc     | Rotation Step | 0     | 5    |

  Scenario Outline: <property> set through the API outside its declared bounds is not taken by the <viewer>
    When user adds a <viewer> viewer
    And user sets "<property>" property of <viewer> viewer to "<set>"
    Then "<property>" property of <viewer> viewer should be "<kept>"
    And no errors should have been logged

    @known-failure:146 @realizes:core#1200
    Examples: The value is written into the property
      | viewer       | property         | set | kept |
      | tree         | initialTreeDepth | 6   | 3    |
      | tree         | fontSize         | 31  | 12   |

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

    @known-failure:166 @ticket-pending
    Examples: The word cloud counts an empty value as a word with no text (`word-cloud-viewer.ts`)
      | viewer     | property       | shown      | reading    |
      | word cloud | wordColumnName | rows shown | word names |

  Scenario Outline: A <viewer> over a column of more than 500 categories keeps its message under a filter that passes fewer than 500 rows
    When user adds a <viewer> viewer with:
      | hierarchyColumnNames | USUBJID |
    Then "<message>" text in <viewer> viewer should be visible
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And "<message>" text in <viewer> viewer should be visible
    And no errors should have been logged

    @known-failure:179 @ticket-pending
    Examples: The category limit is checked only while at least 500 rows pass the filter; under it the viewer draws one branch per USUBJID
      | viewer   | message                                                                                          |
      | tree     | The Tree viewer requires at least one categorical column with fewer than 500 unique categories     |
      | sunburst | The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories |
