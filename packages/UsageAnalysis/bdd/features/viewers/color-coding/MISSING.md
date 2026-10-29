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

## Edge cases (`color-coding-edge-cases.feature`)

Added with `color-coding-edge-cases.feature` (edge-case round of 2026-10-26, localhost). Only the feature was written: no library, binding or code change. Expected behaviour that no help page, TestTrack md, code contract or ticket defines is listed here instead of being written as a scenario.

### Still open

1. **Decided** (2026-10-28):
   - A linked coloring whose source column is removed paints nothing: claimed in
     `color-coding-edge-cases` (core#1139, fixed). The column's menu still shows it as Linked.
   - Bin... takes a number of bins in a range (for example 1 to 100): core#1171, fixed (OK stays
     disabled for 0), claimed in `color-coding-edge-cases`. 10000 bins hangs the tab, so no scenario
     runs it.
   - A coloring picked up from an uncoloured column and applied to a coloured one: the present
     behaviour is right, the coloring is applied as if the source had it on.
2. **Steps requested**: the stop values typed in Edit scheme... (`abc`, stops out of order) have no
   invalid state to read; the conditional editor's "Invalid pattern" mark is the model.
