@guide
Feature: Score compounds with a VlaaiVis summary column
  A guide: the answer to "how do I add a VlaaiVis column to the grid and set it up?". VlaaiVis is a
  summary column that draws each row as a pie: every column it scores is a wedge whose length is the
  value's 0–1 desirability, and columns are grouped into colored sectors. A right-click on any cell,
  Add > Summary Columns > VlaaiVis, appends the column at the end of the grid; dragging its header
  onto the arrow the grid shows at its left edge while a header is dragged makes it the first
  column. A click on the header shows its settings in the Context Panel, under Renderer: Columns
  picks the columns it scores, a sector's header renames or deletes it (its buttons show under the
  pointer), its color square opens a color picker, + Add sector adds one, and a column name dragged
  onto a sector moves it there (Unassigned columns are not drawn). A column's gear opens its score
  function (freeform, gaussian or sigmoid), the number beside it is its weight, and Target min and
  Target max are the shaded ring drawn behind the wedges. A click on a VlaaiVis cell lists the values
  behind it. It starts with a sector for each of the first three numeric columns, which on spgi-100
  are the ID-like CAST Idea ID, Chemical Space X and Chemical Space Y. Demo: spgi-100.

  Scenario: Add a VlaaiVis column, group its columns into sectors and tune their scores
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the context panel is open
    And user opens spgi dataset
    # caption: Right-click any cell and pick Add > Summary Columns > VlaaiVis
    When user picks "Add > Summary Columns > VlaaiVis" from the context menu of the "cell 2 of Id" area of grid
    # caption: The VlaaiVis column is added after the last column
    Then grid should have a "header VlaaiVis" area
    # caption: Drag its header onto the arrow that shows left of the grid: it becomes the first column
    When user drags the "header VlaaiVis" area of grid to the first position
    # caption: Click the VlaaiVis header: the Context Panel shows its settings
    And user clicks on the "header VlaaiVis" area of grid
    # caption: It starts with a sector for each of the first three numeric columns
    Then "CAST Idea ID" VlaaiVis sector should be visible
    # caption: Click Columns to choose the columns it scores
    When user clicks on Columns input in context panel
    # caption: Uncheck CAST Idea ID
    And user toggles the "CAST Idea ID" column in the column list of "Select columns..." dialog
    # caption: Uncheck Chemical Space X
    And user toggles the "Chemical Space X" column in the column list of "Select columns..." dialog
    # caption: Uncheck Chemical Space Y
    And user toggles the "Chemical Space Y" column in the column list of "Select columns..." dialog
    # caption: Uncheck Board ID
    And user toggles the "Board ID" column in the column list of "Select columns..." dialog
    # caption: Check Competition assay
    And user toggles the "Competition assay" column in the column list of "Select columns..." dialog
    # caption: Click OK
    And user clicks on OK button in "Select columns..." dialog
    # caption: Seven columns are scored now; the ones no sector holds are listed under Unassigned columns
    Then Columns input in context panel should contain text "(7) Average Mass, TPSA"
    # caption: Point at the CAST Idea ID sector to show its buttons
    When user hovers over header of "CAST Idea ID" VlaaiVis sector
    # caption: Click the pencil to rename the sector
    And user clicks on "Rename sector" icon in "CAST Idea ID" VlaaiVis sector
    # caption: Name it Physchem
    And user types "Physchem" into Name input in "Rename sector" dialog
    # caption: Click OK
    And user clicks on OK button in "Rename sector" dialog
    # caption: Drag TPSA from Unassigned columns onto Physchem
    And user drags "TPSA" text in "Unassigned columns" VlaaiVis sector to "Physchem" VlaaiVis sector
    # caption: Drag NIBR logP onto Physchem too
    And user drags "NIBR logP" text in "Unassigned columns" VlaaiVis sector to "Physchem" VlaaiVis sector
    # caption: Both are wedges of the Physchem sector now
    Then "NIBR logP" VlaaiVis column in "Physchem" VlaaiVis sector should be visible
    # caption: Point at the Chemical Space X sector
    When user hovers over header of "Chemical Space X" VlaaiVis sector
    # caption: Click the pencil
    And user clicks on "Rename sector" icon in "Chemical Space X" VlaaiVis sector
    # caption: Name it Potency
    And user types "Potency" into Name input in "Rename sector" dialog
    # caption: Click OK
    And user clicks on OK button in "Rename sector" dialog
    # caption: Drag Competition assay onto Potency
    And user drags "Competition assay" text in "Unassigned columns" VlaaiVis sector to "Potency" VlaaiVis sector
    # caption: Competition assay is Potency's wedge now
    Then "Competition assay" VlaaiVis column in "Potency" VlaaiVis sector should be visible
    # caption: Point at the empty Chemical Space Y sector
    When user hovers over header of "Chemical Space Y" VlaaiVis sector
    # caption: Click the trash can to delete it
    And user clicks on "Delete sector" icon in "Chemical Space Y" VlaaiVis sector
    # caption: The sector is gone
    Then "Chemical Space Y" VlaaiVis sector should be absent
    # caption: Click + Add sector for a new one
    When user clicks on "+ Add sector" link
    # caption: Name it Size
    And user types "Size" into Name input in "New sector" dialog
    # caption: Click OK
    And user clicks on OK button in "New sector" dialog
    # caption: Drag Num Heavy Atoms onto Size
    And user drags "Num Heavy Atoms" text in "Unassigned columns" VlaaiVis sector to "Size" VlaaiVis sector
    # caption: Size has its first column
    Then "Num Heavy Atoms" VlaaiVis column in "Size" VlaaiVis sector should be visible
    # caption: Click Potency's color square to change its color
    When user clicks on color of "Potency" VlaaiVis sector
    # caption: Pick red
    And user picks the color "#d62728" in the color picker
    # caption: Click OK: Potency's wedges turn red
    And user clicks on OK button
    # caption: Click TPSA's gear to edit its score function
    And user clicks on "Settings" icon in "TPSA" VlaaiVis column
    # caption: Switch the curve to gaussian
    And user selects "gaussian" in Mode input in "TPSA" dialog
    # caption: Make 75 the most desirable TPSA
    And user types "75" into Mean input in "TPSA" dialog
    # caption: Make the bell narrower
    And user types "12" into Sigma input in "TPSA" dialog
    # caption: Click OK
    And user clicks on OK button in "TPSA" dialog
    Then the "TPSA" dialog should close
    # caption: Halve NIBR logP's weight: its wedge gets narrower
    When user types "0.5" into "NIBR logP" VlaaiVis column
    # caption: Set Target min, the inner edge of the shaded ring
    And user types "0.6" into "Target min" input in context panel
    # caption: and Target max, its outer edge
    And user types "0.8" into "Target max" input in context panel
    # caption: Click a VlaaiVis cell: the Context Panel lists the values behind it
    And user clicks on the "cell 3 of VlaaiVis" area of grid
    # caption: The values are grouped by sector, with the unassigned columns last
    Then context panel should contain text "Physchem"
