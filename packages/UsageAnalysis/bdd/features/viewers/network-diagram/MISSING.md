# Network diagram — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against
`files/TestTrack/Viewers/NetworkDiagram/network-diagram.md` and the viewer's "Automation surface"
(`core/client/d4/lib/src/viewers/network_diagram/CLAUDE.md`). The status reports nodes (areas, names,
`rows of node`) but nothing about edges, the empty canvas or the view's scale.

## Still open

- **Clicking an edge selects the rows behind it** (md "Selection gestures", step 6). No `edge "<from>
  -> <to>"` area and no `rows of edge` reading. Needs both, from vis's edge geometry
  (`Network.getConnectedNodes` / the edge's `from`/`to` positions through `canvasToDOM`, the middle
  of the edge as the pointing target), and the rows the click selects (the same rows the edge
  aggregates):

  ```gherkin
  Scenario: Clicking an edge selects only the rows behind that edge
    Given user clears the row selection
    When user clicks on the "edge \"F -> false\"" area of network diagram viewer
    Then the "rows of edge \"F -> false\"" reading of network diagram viewer should be <measured>
    And <measured> rows should be selected
    And no errors should have been logged
  ```

  Select Edges On Click off (md step 2 of "Misc") then needs the same area to claim the click selects
  nothing.

- **Double-click on empty canvas** (md step 7: the first click clears the selection, the second zooms
  to fit). No `empty space` area (a point no node or edge covers) and no `view scale` reading
  (`Network.getScale()`), so neither half can be claimed:

  ```gherkin
  Scenario: A double-click on empty canvas clears the selection and fits the view
    Given user clicks on the "node \"F\"" area of network diagram viewer
    And user scrolls the mouse wheel up over the "view" area of network diagram viewer
    And user remembers the "view scale" reading of network diagram viewer
    When user double-clicks on the "empty space" area of network diagram viewer
    Then no rows should be selected
    And the "view scale" reading of network diagram viewer should not be as remembered
    And no errors should have been logged
  ```

- **Arrow heads.** `network-diagram.feature` claims only the `arrows` setting and a repaint (its own
  description says so). Needs an `arrow heads` reading — how many edges vis draws with an arrow
  (`edges.get()` whose `arrows` option is set) — so Show Arrows = to can be claimed as 4 arrow heads
  on demog-1000 and 0 after it is switched off.
