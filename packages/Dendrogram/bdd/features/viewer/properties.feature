@journey @realizes:dendrogram.cp.viewer-from-newick-prop
Feature: The Dendrogram viewer's properties from the property panel
  A four-leaf table (leaf a to d, value 1 to 4) tagged with two newick trees — .newick balanced,
  .newick-alt left-leaning — and a Dendrogram viewer given the balanced tree in its Newick property
  and the leaf column as its node column. Its property panel shows the Data, Style and Behavior
  categories. Newick wins over Newick Tag, and with Newick cleared the tag it names is drawn; the node
  column binds the grid's current row to a node; Color and Color Aggr Type reach the main style; Line
  Width, Node Size and Show Grid reach the main style; the colors reach the main, light, current and
  selection styles. Show Tooltip is not claimed.

  Mouse Over Color is not applied to its style, and Show Labels, Font, Step and Step Zoom are held by
  neither the renderer nor the placer, so the tree draws the same whatever they are set to. The last
  five scenarios state what each should reach and are known failures; each tag goes when that
  property reaches the tree.

  Background:
    Given user is logged in
    And user opens a table "leaves4" with:
      | leaf | value |
      | a    | 1     |
      | b    | 2     |
      | c    | 3     |
      | d    | 4     |
    And user sets tag ".newick" of the table to "((a:1,b:1):1,(c:1,d:1):1);"
    And user sets tag ".newick-alt" of the table to "(((a:1,b:1):1,c:1):1,d:1);"
    And user adds Dendrogram viewer with:
      | newick         | ((a:1,b:1):1,(c:1,d:1):1); |
      | nodeColumnName | leaf                       |
    Then the "newick" reading of Dendrogram viewer should be "((a:1,b:1):1,(c:1,d:1):1);"
    And the "leaves" reading of Dendrogram viewer should be "a, b, c, d"
    When user picks "Properties..." from the context menu of Dendrogram viewer

  Scenario: The property panel shows the viewer's categories
    Then Data category should be visible
    And Style category should be visible
    And Behavior category should be visible
    And "Newick Tag" property should be visible
    When user expands Style category
    And user expands Behavior category
    Then "Line Width" property should be visible
    And "Show Tooltip" property should be visible
    And no errors should have been logged

  Scenario: Newick wins over Newick Tag; with Newick cleared the named tag is drawn
    When user selects ".newick-alt" in "Newick Tag" property
    Then the "newick" reading of Dendrogram viewer should be "((a:1,b:1):1,(c:1,d:1):1);"
    When user clears Newick property
    And user presses Enter in Newick property
    Then the "newick" reading of Dendrogram viewer should be "(((a:1,b:1):1,c:1):1,d:1);"
    And the "leaves" reading of Dendrogram viewer should be "a, b, c, d"
    When user selects "" in "Newick Tag" property
    Then the "newick" reading of Dendrogram viewer should be "((a:1,b:1):1,(c:1,d:1):1);"
    And no errors should have been logged

  Scenario: The node column binds the grid's current row to a node
    When user clicks on the "cell 2 of leaf" area of grid
    Then the "current node" reading of Dendrogram viewer should be "b"
    When user selects "value" in Node property
    And user clicks on the "cell 3 of leaf" area of grid
    Then the "current node" reading of Dendrogram viewer should be ""
    When user selects "leaf" in Node property
    And user clicks on the "cell 4 of leaf" area of grid
    Then the "current node" reading of Dendrogram viewer should be "d"
    And no errors should have been logged

  Scenario: Color and Color Aggr Type color the tree by each aggregation
    When user selects "value" in Color property
    And user selects "avg" in "Color Aggr Type" property
    Then the "color coding" reading of Dendrogram viewer should be "avg of value"
    When user selects "min" in "Color Aggr Type" property
    Then the "color coding" reading of Dendrogram viewer should be "min of value"
    When user selects "max" in "Color Aggr Type" property
    Then the "color coding" reading of Dendrogram viewer should be "max of value"
    When user selects "med" in "Color Aggr Type" property
    Then the "color coding" reading of Dendrogram viewer should be "med of value"
    When user selects "count" in "Color Aggr Type" property
    Then the "color coding" reading of Dendrogram viewer should be "count of value"
    And no errors should have been logged

  Scenario: Line Width, Node Size and Show Grid reach the main style
    When user enters "0" into "Line Width" property
    Then the "line width" reading of Dendrogram viewer should be 0
    When user enters "16" into "Line Width" property
    Then the "line width" reading of Dendrogram viewer should be 16
    When user enters "2.5" into "Line Width" property
    Then the "line width" reading of Dendrogram viewer should be 2.5
    When user enters "0" into "Node Size" property
    Then the "node size" reading of Dendrogram viewer should be 0
    When user enters "16" into "Node Size" property
    Then the "node size" reading of Dendrogram viewer should be 16
    When user enters "4" into "Node Size" property
    Then the "node size" reading of Dendrogram viewer should be 4
    When user checks "Show Grid" property
    Then the "show grid" reading of Dendrogram viewer should be "true"
    When user unchecks "Show Grid" property
    Then the "show grid" reading of Dendrogram viewer should be "false"
    And no errors should have been logged

  Scenario: The style colors reach their styles
    When user clicks on value of "Main Color" property
    And user picks the color "#d62728" in the color picker dialog
    And user clicks on label of "Step Zoom" property
    Then the "main color" reading of Dendrogram viewer should be "#d62728"
    When user clicks on value of "Light Color" property
    And user picks the color "#d62728" in the color picker dialog
    And user clicks on label of "Step Zoom" property
    Then the "light color" reading of Dendrogram viewer should be "#d62728"
    When user clicks on value of "Current Color" property
    And user picks the color "#d62728" in the color picker dialog
    And user clicks on label of "Step Zoom" property
    Then the "current color" reading of Dendrogram viewer should be "#d62728"
    When user clicks on value of "Selections Color" property
    And user picks the color "#d62728" in the color picker dialog
    And user clicks on label of "Step Zoom" property
    Then the "selections color" reading of Dendrogram viewer should be "#d62728"
    And no errors should have been logged

  @known-failure
  Scenario: Mouse Over Color reaches its style
    When user clicks on value of "Mouse Over Color" property
    And user picks the color "#d62728" in the color picker dialog
    And user clicks on label of "Step Zoom" property
    Then the "mouse over color" reading of Dendrogram viewer should be "#d62728"

  @known-failure
  Scenario: Show Labels makes the tree draw its labels
    When user checks "Show Labels" property
    Then the "labels drawn" reading of Dendrogram viewer should be "true"

  @known-failure
  Scenario: Font reaches the labels
    When user enters "12pt monospace" into Font property
    Then the "label font" reading of Dendrogram viewer should be "12pt monospace"

  @known-failure
  Scenario: Step sets the leaf row spacing
    When user enters "40" into Step property
    Then the "row step" reading of Dendrogram viewer should be "40"

  @known-failure
  Scenario: Step Zoom sets the zoom step
    When user enters "2" into "Step Zoom" property
    Then the "zoom step" reading of Dendrogram viewer should be "2"
