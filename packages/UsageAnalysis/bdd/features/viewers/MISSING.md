# Cross-viewer features — what cannot be said yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026 for the features that sit in
this folder rather than in a viewer's own: `viewers-docking.feature`,
`viewers-layout-download.feature` (added in that round: View > Layout > Download gives a `.layout`
file holding the docked viewers) and `molecule-categories.feature`.

## Still open

- **The splitter between two docked viewers** (`viewers-docking-ui.md`: "drag the splitter between
  them"). dock-spawn's splitter bars carry no name and report no box, so neither a drag nor "the
  left viewer got wider" can be aimed. Needs a name per splitter (`dock-splitter-<n>`, or named by
  the two panels it divides) in `core/client/libs/dock_spawn`; the claim is then
  `histogram viewer should be wider than before` from a `user drags {element} by N and 0 pixels`.
- **Dropping a downloaded `.layout` back onto the view** (`viewers-docking-ui.md`, last step): a
  file dragged in from the desktop. `user drops the {string} file of the project onto {element}`
  takes a fixture file, not the file just downloaded; a step that drops the downloaded file would
  close the round trip.
- **The range sliders over molecule categories on the box plot and the trellis plot**
  (`rendering-structures-on-the-axes-ui.md`, #2330). The bar chart and the scatter plot are walked
  in `molecule-categories.feature`. The box plot reports no slider area (`box_plot/CLAUDE.md`), so
  its category slider cannot be aimed at. The trellis plot reports its `x scroll handle`, and a drag
  of it moves the categories on screen, but every reading stays the same (`x categories shown` 5,
  `x scroll handle share` 0.05 on spgi-100 by Structure) and the label areas are named by the
  molblock text — no phrase can spell them. Needs a `first x category` reading (its row index or
  position) on the trellis, and `x range slider` / handle areas on the box plot:

  ```gherkin
  Scenario: Dragging the trellis plot's category scroll handle over the molecules moves the window
    When user adds a trellis plot viewer with:
      | xColumnNames | Structure |
    Then the "x scroll handle share" reading of trellis plot viewer should be 0.05
    When user drags the "x scroll handle" area of trellis plot viewer by 100 pixels to the right
    Then the "first x category" reading of trellis plot viewer should be higher than before
    And the "blank cells" reading of trellis plot viewer should be 0
    And trellis plot viewer should report no error
  ```
