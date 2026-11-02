# Tile viewer — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against
`files/TestTrack/Viewers/TileViewer/tile-viewer-ui.md` and the viewer's "Automation surface"
(`core/client/d4/lib/src/viewers/tile_viewer/CLAUDE.md`). New this round:
`tile-viewer-designed-rebind.feature` (a card designed on spgi-100, then Table → demog-1000) and
`tile-viewer-lane-scroll.feature` (GROK-20912, `@known-failure`).

## Still open

- **Colour coding on tiles** (`tile-viewer-ui.md`, "Colour coding on tiles"). A card field takes the
  grid column's colour as the input's inline `backgroundColor` (`sketch_field.dart:63-75`,
  `input.setColor`), and is cleared back when the column loses its colouring. Nothing reads it: the card
  is DOM, so the area colour steps (`the "field AGE of row 1" area ... should contain the color`)
  fail with "has no canvas", and no element phrase reaches a card's input (`"AGE" input in tile viewer`
  resolves to nothing). Needs a `background of <COL> of row <r>` reading (computed background as
  `#rrggbb`, `''` when transparent — the forms viewer's `background of <COL> of <label>` is the model):

  ```gherkin
  Scenario: A linear AGE and a categorical RACE colour the fields, and switching them off clears them
    When user colors "AGE" column linearly from "#0000ff" to "#ff0000"
    Then the "background of AGE of row 1" and "background of AGE of row 3" readings of tile viewer should differ
    When user colors "RACE" column categorically:
      | Caucasian | #1f77b4 |
      | Other     | #ff7f0e |
    Then the "background of RACE of row 1" reading of tile viewer should be "#1f77b4"
    And the "background of RACE of row 2" reading of tile viewer should be "#ff7f0e"
    When user removes the coloring of "AGE" column
    And user removes the coloring of "RACE" column
    Then the "background of AGE of row 1" reading of tile viewer should be ""
    And the "background of RACE of row 1" reading of tile viewer should be ""
    And no errors should have been logged
  ```

- **Picking the lanes column in the property panel** (`tile-viewer-ui.md`). The Lanes column picker
  is a column combo in the property grid; every lanes scenario sets the property through the API.
  Needs a library step that opens a column-typed property row's picker and chooses a column
  (`user picks {string} for the {string} property of {widget}`), then:

  ```gherkin
  Scenario: Choosing RACE in the Lanes picker regroups the cards
    When user picks "RACE" for the "Lanes Column Name" property of tile viewer
    Then the "lanes" reading of tile viewer should be 4
    And the "lane names" reading of tile viewer should be "Asian, Black, Caucasian, Other"
    And no errors should have been logged
  ```

## Edge cases (now in `tile-viewer.feature`, `tile-viewer-form-editor.feature` and `tile-viewer-known-failures.feature`)

Written for the guard suite of edge cases of 26 Oct 2026 (its checks are now in
`tile-viewer.feature`, and `@edge-cases` scenarios of their own after the master scenarios of
`tile-viewer-form-editor.feature`, and the known failure in `tile-viewer-known-failures.feature`): data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

None.

#### 2. Decided

1. **The lanes column removed** is dropped from the settings, as clearing it does, and no other
   column is picked (`TileViewerLook.auto` keeps only a lanes name the table has). This check, with
   no rows and a renamed lanes column, is requested from ApiTests in `../MISSING.md`.
2. **A lanes column with empty values** has an `(empty)` lane: core#1217, still open, requested
   from ApiTests in `../MISSING.md`. That lane offered as `(empty)` is written in
   `tile-viewer-known-failures.feature` (back into `tile-viewer-lanes.feature` once fixed), core#1217, a
   known failure (the lane has no name today, so there is no `(empty)` lane to drop on);
   the drop and the empty value it writes are claimed once the lane exists.
3. **Escape in the form designer** after deleting a field does nothing: the designer is a view, not
   a dialog. Right as it is (2026-10-30), claimed in `tile-viewer-form-editor.feature`. Escape in the
   Edit dialog (the EDIT icon of the ribbon) cancels it: not written.

#### 3. Questions for owners

None.

**Decided (2026-10-30):** a lanes column with more than 20 categories cannot be picked as Lanes; closing
the form designer's tab after deleting a field works as it does today. Neither has a scenario.

**Decided:** a drag between the lanes of a calculated column writes the lane's value: core#1237 was
closed as won't fix (2026-10-30), not tested.

#### 3. API or Dart tests to write

1. Dart regTest: the Data Manipulation loop (`viewers_test.dart:904-950`) runs for
   `ViewerDescriptor.all.take(30)` and skips the tile viewer (index 30); drop the `take(30)`.
2. ApiTests: Row Source Selected with nothing selected shows 0 rows and 0 tiles; selecting the
   Asian rows brings 15 rows, with row 10 in lane Asian; nothing logged.
3. ApiTests: a one-row table with lanes on its group column shows 1 row and 1 tile, in lane g1.

Moved out of BDD (no user gesture); requested from ApiTests.

### Edge cases: Won't do

- core#1088 (a designer edit reaches one lane only): already pinned in
  `tile-viewer-persistence.feature:71`.

## Closed issues the regression round of 26 Oct 2026 could not claim

- **A value the card cuts off shows whole in a tooltip (public#1593, the value half)** regressed: no tooltip
  shows; filed as core#1151. The caption half is claimed in `tile-viewer-mirroring.feature`.
