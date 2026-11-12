@journey @realizes:dendrogram.viewer.newick-property
Feature: The Dendrogram viewer drawn from a newick it is given
  The Dendrogram viewer draws the tree of a newick: its own Newick property first, else the table tag
  its Newick Tag property names, else the table's .newick tag. A viewer added with a literal newick
  draws that tree whatever the table's tags hold; cleared, it falls back to the .newick tag; and the
  Newick Tag editor of the context panel offers the table's dot-tags, attached when the viewer joined
  the table, so picking .newick-alt draws that tree. Its style properties restyle the tree, and a
  color column colors the leaves by the column's own color coding. Translated from TestTrack
  Dendrogram/dendrogram-viewer-from-newick-prop.md and the package's
  playwright/dendrogram-viewer-from-newick-prop.test.ts.

  The trees are told apart by the viewer's own readings: the newick it draws and its leaves in
  drawing order. The table's .newick tag is the mirror of the given tree, so a viewer that drew the
  tag instead of its property would read "d, c, b, a". What the style properties change is claimed
  as a repaint of the tree's canvas (the shapes are the renderer's, tested in the package). A tree
  whose root has no branch length, as most newick strings have it, drew nothing at all until
  2026-10-13 (markup.ts summed an undefined length into the tree's height).

  Show Labels, Step, Step Zoom, Font, Show Tooltip, Mouse-over Color and Color Aggr Type are
  declared but nothing in the viewer reads them (dendrogram.ts, 2026-10-13): Show Labels is claimed
  as a known failure — switched on, the tree draws no leaf names and does not repaint at all.

  Not translated: the node column's binding of the current, mouse-over and selected rows — the
  viewer reports no hit area for its nodes, so no gesture can aim at one.

  Background:
    Given user is logged in
    And user opens a table "leaves4" with:
      | leaf | value |
      | a    | 1     |
      | b    | 2     |
      | c    | 3     |
      | d    | 4     |
    # the table view makes a row current a second after it opens, and the tree then marks that row's leaf
    And user makes row 1 current
    And user sets the ".newick" tag of the table to "((d:1,c:1):1,(b:1,a:1):1);"
    And user sets the ".newick-alt" tag of the table to "(((a:1,b:1):1,c:1):1,d:1);"
    And user adds a Dendrogram viewer with:
      | newick         | ((a:1,b:1):1,(c:1,d:1):1); |
      | nodeColumnName | leaf                       |

  Scenario: The newick the viewer is given wins over the table's tags
    Then Dendrogram viewer should be visible
    And the "newick" reading of Dendrogram viewer should be "((a:1,b:1):1,(c:1,d:1):1);"
    And the "leaves" reading of Dendrogram viewer should be "a, b, c, d"
    And Dendrogram viewer should be painted
    And no errors should have been logged

  @known-failure:52
  Scenario: Show Labels draws the leaf names
    When user sets "showLabels" property of Dendrogram viewer to "true"
    Then Dendrogram viewer should have more ink than before
    When user sets "showLabels" property of Dendrogram viewer to "false"

  Scenario: Without its own newick the viewer draws the tag Newick Tag names, else .newick
    When user sets "newick" property of Dendrogram viewer to ""
    Then the "newick" reading of Dendrogram viewer should be "((d:1,c:1):1,(b:1,a:1):1);"
    And the "leaves" reading of Dendrogram viewer should be "d, c, b, a"
    When user clicks on settings icon of Dendrogram viewer
    Then context panel should be visible
    When user selects ".newick-alt" in "Newick Tag" property in context panel
    Then the "newick" reading of Dendrogram viewer should be "(((a:1,b:1):1,c:1):1,d:1);"
    And the "leaves" reading of Dendrogram viewer should be "a, b, c, d"
    When user sets "newickTag" property of Dendrogram viewer to ""
    And user sets "newick" property of Dendrogram viewer to "((a:1,b:1):1,(c:1,d:1):1);"
    Then the "newick" reading of Dendrogram viewer should be "((a:1,b:1):1,(c:1,d:1):1);"
    And no errors should have been logged

  Scenario: A color column colors the leaves by its color coding
    When user colors "value" column linearly from "#0000ff" to "#ff0000"
    And user sets "colorColumnName" property of Dendrogram viewer to "value"
    Then Dendrogram viewer should have repainted
    And Dendrogram viewer should be painted in at least 3 colors
    When user sets "colorColumnName" property of Dendrogram viewer to ""
    Then Dendrogram viewer should have repainted
    And no errors should have been logged

  Scenario: The style properties restyle the tree
    When user sets "lineWidth" property of Dendrogram viewer to "8"
    Then Dendrogram viewer should have repainted
    When user sets "nodeSize" property of Dendrogram viewer to "12"
    Then Dendrogram viewer should have repainted
    When user sets "showGrid" property of Dendrogram viewer to "true"
    Then Dendrogram viewer should have repainted
    When user sets "mainColor" property of Dendrogram viewer to "#ff0000"
    Then Dendrogram viewer should have repainted
    And no errors should have been logged
