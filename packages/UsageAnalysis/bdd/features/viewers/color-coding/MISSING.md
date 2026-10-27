# Color coding editor — what the features cannot say yet

Written during the color scheme editor round of 14 Oct 2026 (TestTrack `Viewers/color-coding-ui.md`
steps 9-12, the help's "Color code columns", and the closed Jira and GitHub tickets on color coding),
reviewed on 18 Oct 2026 against the library's rules (`public/libraries/bdd/CLAUDE.md`).

## Resolved

- **A category's colour in the categorical editor (public#3413).** The `item` kind has a `color bar` part
  (`.d4-legend-item-bar`, the bar that opens the picker), so `color-coding-editor-conditional-and-categorical`
  picks a colour for "M" with the existing picker step and claims the column's categorical colour and the
  colour the grid paints the cell, before Off clears the column.
- **The text flag unset.** Switching Apply to back to background stored the tag as the string "null"
  (`GridColumn.isTextColorCoded` set it without `removeIfNull`, `grid_column.dart`); it is now removed, pinned
  by ApiTests `Grid: Color Coding` "colorCoding.text", and `the text of {string} column should not be
  color-coded` claims both directions of Apply to in `color-coding-editor-linear`.
- **A saved layout bringing the text coloring back (GROK-13668).** What a layout stores: ApiTests `Grid: Color
  Coding` "colorCoding.api" now sets the text flag on SEX and claims it after the layout is applied to a new
  view; the scenario is gone from the feature.
- **The conditional rows and the stops of "Edit color scheme".** A map input's row editors are named ("Key N" /
  "Value N", `map_input.dart`); the conditional editor calls its rows "Pattern N" with a "Color N" bar, the scheme
  editor its stops "Stop N" with a "Stop N color" bar, and a colour input's box is a button named by its caption
  ("Below min"), all in core d4. The `text input` kind takes a map input's named row editor. The features type a
  pattern that takes cells over, a single value that colours the cells equal to it (public#1967) and an invalid
  pattern with its "Invalid pattern" tooltip, pick a stop's colour and Below min's.

## Still open

Nothing planned.

## Won't do

- **Only one outlier picker open (GROK-17738)**: the Below min picker opens over the Above max bar, so a person
  cannot open the second without closing the first; there is no state to claim.
- **The text colour of a cell** (public#2784): what the renderer paints is the renderer's; the features claim
  the column's text flag and the colour the grid resolved.
- **No Cancel in the editors**: intended (Olesia Pavlenko, 14 Oct 2026); `color-coding-ui.md` step 10 and
  GROK-8444 describe a dialog that no longer exists. **Off keeping the text flag**: intended too.
- **What a layout or a project stores for a coloring**: ApiTests `Grid: Color Coding`.
- **Pixel checks of a renderer, a correlation plot's mirrored cells (GROK-20919)**: renderer.
- **The legend's own picker, its Cancel and a second pick (github-3132)**: `legend/legend-across-viewers.feature`;
  a conditional coloring's legend: `guides/scatter-plot-conditional-colors.feature`.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed colour coding issues into `color-coding-regressions.feature`.

### Regressions found 26 Oct 2026

- **A colour-coded category column paints the bar chart's category names (public#3621).** Measured twice
  on localhost: on demog-1000 with RACE coloured categorically (Caucasian #FF0000, Other #00FF00) and
  a bar chart split by RACE, the `y axis` area (the category names) is painted #FF0000 — 7009 px of it
  in the 75×200 area, every one of the four names on a red background, the green "Other" included,
  while the grid cells are coloured as set. `_renderLabel` in `bar_chart_core.dart` sets the style's
  `backColor` to transparent, but the category renderer it calls still paints the cell. Expected (the
  issue): no colouring on the category names.

  ```gherkin
  @realizes:public#3621
  Scenario: A colour-coded category column does not paint the bar chart's category names (public#3621)
    Given user is logged in
    And user opens demog-1000 dataset
    When user colors "RACE" column categorically:
      | Caucasian | #FF0000 |
      | Other     | #00FF00 |
    Then the "cell 1 of RACE" area of grid should contain the color "#FF0000"
    When user adds a bar chart viewer with:
      | Split | RACE |
    Then bar chart viewer should have a "y axis" area
    And the "y axis" area of bar chart viewer should not contain the color "#FF0000"
    And the "y axis" area of bar chart viewer should not contain the color "#00FF00"
    And no errors should have been logged
  ```

### Questions for owners

- **Color Coding > Apply Coloring in another table (public#1887's note).** Pick Up / Apply > Pick Up on
  demog-1000's RACE, then Pick Up / Apply > Apply on demog-100's RACE, carries the colouring (claimed).
  Color Coding > Pick Up Coloring on demog-1000's RACE, then the header menu of demog-100's RACE: Apply
  Coloring is disabled with "Pick up coloring first" (seen 26 Oct 2026). The issue's note says applying
  the colour coding should work from a style picked up in a different dataframe: is the Color Coding
  pair meant to be per table?
