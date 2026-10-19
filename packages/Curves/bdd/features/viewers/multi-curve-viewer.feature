@journey
Feature: The Multi curve viewer from the Add Viewer dialog, its curve columns, Use in Trellis
  Translated from TestTrack Curves/multi-curve-viewer.md, scenarios 1, 2, 3, 4 (the viewer's part), 5
  and 7. The MultiCurveViewer card of the Add Viewer dialog is enabled on a table with curve columns and
  disabled on one without (GROK-19014, GROK-14069); added, the viewer takes one of the five curve
  columns and draws the curves of the current row, the Curves editor's All takes all five, the selected
  rows add theirs; Merge Cell Series draws one curve per cell (GROK-19028); General > Use in Trellis
  puts it into a trellis as the inner viewer (GROK-19025); a closed viewer leaves no error behind
  (Curves 1.13.1). The curves table is ready once the Context Panel shows the value of its current
  cell: until then the platform may still replace the viewer's properties opened from its menu.

  What the viewer draws is its "curves" reading. Row 1 of the curves table holds 2 + 2 + 1 + 1 + 1
  curves in its five columns. That merging leaves the grid's own cell alone is the package's
  multi-curve-viewer test; the legend is the guide guides/curves-legend.feature. With all five columns, every
  point lies within the axes the viewer draws (GROK-20661: a range one cell declares for its own chart no
  longer narrows the overlay).

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed

  Scenario: On a table without curves the MultiCurveViewer card is disabled
    Given user opens demog-1000 dataset
    Then "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "curve" into viewer gallery search in "Add Viewer" dialog
    Then first "MultiCurveViewer" card in "Add Viewer" dialog should be disabled
    When user presses Escape
    Then "Add Viewer" dialog should be absent
    And no errors should have been logged

  Scenario: On the curves table the card adds the viewer with one curve column, and All takes five
    Given user closes all views
    And user opens curves dataset
    Then the context panel should show the current cell
    And "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "curve" into viewer gallery search in "Add Viewer" dialog
    Then first "MultiCurveViewer" card in "Add Viewer" dialog should be enabled
    When user clicks on first "MultiCurveViewer" card in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And the open tableview should have 1 MultiCurveViewer viewer
    And the "curves" reading of MultiCurveViewer viewer should be 2
    When user picks "Properties..." from the viewer menu of MultiCurveViewer viewer
    Then "Curves" property in context panel should contain text "1 / 5"
    When user clicks on "..." button in "Curves" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    Then "5 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    And user picks "Properties..." from the viewer menu of MultiCurveViewer viewer
    Then "Curves" property in context panel should contain text "5 / 5"
    And the "curves" reading of MultiCurveViewer viewer should be 7
    And the "points outside the axes" reading of MultiCurveViewer viewer should be 0
    And no errors should have been logged

  Scenario: The selected rows add their curves, Merge Cell Series draws one curve per cell
    When user selects the first 3 rows
    Then the "curves" reading of MultiCurveViewer viewer should be 21
    When user expands "Misc" category in context panel
    And user checks "Merge Cell Series" property in context panel
    Then the "curves" reading of MultiCurveViewer viewer should be 15
    When user unchecks "Merge Cell Series" property in context panel
    Then the "curves" reading of MultiCurveViewer viewer should be 21
    When user clears the row selection
    Then the "curves" reading of MultiCurveViewer viewer should be 7
    And no errors should have been logged

  Scenario: A closed viewer leaves no error when rows are clicked afterwards
    When user hovers over MultiCurveViewer viewer
    And user clicks on close icon of MultiCurveViewer viewer
    Then the open tableview should have 0 MultiCurveViewer viewers
    When user selects the first 3 rows
    And user makes row 4 current
    Then no errors should have been logged

  Scenario: Use in Trellis puts the viewer into a trellis plot as its inner viewer
    Given user adds MultiCurveViewer viewer
    When user picks "General | Use in Trellis" from the viewer menu of MultiCurveViewer viewer
    Then the open tableview should have 1 trellis plot viewer
    And the "inner viewer type" reading of trellis plot viewer should be "MultiCurveViewer"
    And trellis plot viewer should report no error
    And no errors should have been logged
