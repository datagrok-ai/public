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

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed tile viewer issues into `tile-viewer-regressions.feature`.

### Regressions found 26 Oct 2026

- **A value the card cuts off shows whole in a tooltip (public#1593, the value half).** Measured twice on
  localhost (a BDD run and a headless probe with a real 1.5 s hover): on a two-row table whose `note`
  is 120 characters long, the card's field shows "The first value of" and hovering `field note of row 1`
  shows no tooltip (no visible `.d4-tooltip`). The caption half of the same fix works (claimed in
  `tile-viewer-regressions.feature`). The field's tooltip is set in `sketch_field.dart` `setRow`
  (`tooltip.setOn(input.root, …)`, check: text width > `input.defaultSize.width`, 100 px), moved there
  from `sketch_form.dart` by 9e602c84f3 (Sep 2024). Expected (the issue): the whole value in a tooltip.

  ```gherkin
  @realizes:public#1593
  Scenario: A value the card cuts off shows whole in a tooltip (public#1593)
    Given user is logged in
    And user opens a table "long-values" with:
      | id | note                                                                                                                     |
      | 1  | The first value of this column is long enough that no field of a card can show it whole, so its end is cut off here. |
      | 2  | short                                                                                                                    |
    And user adds a tile viewer
    Then the "note of row 1" reading of tile viewer should be "The first value of this column is long enough that no field of a card can show it whole, so its end is cut off here."
    When user hovers over the "field note of row 1" area of tile viewer
    Then tooltip should be visible
    And tooltip should contain text "so its end is cut off here."
    When user moves the pointer away from tile viewer
    Then no errors should have been logged
  ```

### Steps requested

- **A calculated column's formula edited in the Formula pane updates the card (public#3344).** The
  Formula pane's editor is the noun `formula pane editor`, registered by PowerPack's bindings, not in
  this project ("neither a registered element nor <qualifier> <kind>" at compile). Needs the noun in the
  library (or this scenario in PowerPack's project), then:

  ```gherkin
  @realizes:public#3344
  Scenario: A calculated column's formula edited in the Formula pane updates the card (public#3344)
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer
    When user removes "DEMOG" column
    And user removes "STARTED" column
    And user adds a calculated column "cc" with formula "${AGE} + 5"
    Then the "cc of row 1" reading of tile viewer should be "31"
    Given the context panel is open
    When user clicks on the "header cc" area of grid
    Then the context panel should show "cc"
    Given Formula pane in context panel is expanded
    Then formula pane editor in Formula pane in context panel should hold the code "${AGE} + 5"
    When user types "${AGE} + 10" into formula pane editor in Formula pane in context panel
    And user clicks on Apply button in Formula pane in context panel
    Then the value of "cc" column in row 1 should be "36"
    And the "cc of row 1" reading of tile viewer should be "36"
    And no errors should have been logged
  ```
  (Not run: the "31" reading was not reached in this project; the numbers are AGE 26 + 5 and + 10.)
