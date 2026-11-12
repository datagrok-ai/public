# Color coding editor — what the features cannot say yet

Written during the color scheme editor round of 14 Oct 2026 (TestTrack `Viewers/color-coding-ui.md`
steps 9-12, the help's "Color code columns", and the closed Jira and GitHub tickets on color coding),
reviewed on 18 Oct 2026 against the library's rules (`public/libraries/bdd/CLAUDE.md`) and in the BDD PR
round of 8 Nov 2026.

## Resolved

- **A category's colour in the categorical editor (public#3413)**: the `color bar` part of an `item`, claimed in
  `color-coding-editor-conditional-and-categorical`.
- **The text flag unset**: ApiTests `Grid: Color Coding` "colorCoding.text"; both directions of Apply to in
  `color-coding-editor-linear`.
- **A saved layout bringing the text coloring back (GROK-13668)**: ApiTests `Grid: Color Coding` "colorCoding.api".
- **The conditional rows and the stops of "Edit color scheme"** are named ("Pattern N", "Stop N", "Below min"):
  claimed in the editor features.

## Still open

- **core#1153** — a colour-coded Split column no longer paints the bar chart's category names (public#3621
  regressed); claim it in `bar-chart/bar-chart.feature` once fixed.
- **Needs a ruling**: a one-value column coloured linearly or conditionally paints nothing.

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
- **Pick Up / Apply Coloring across grids**: per grid by design (Olesia, 27 Oct 2026); across tables it is Pick Up /
  Apply (public#1887).
- **Bin... with 0 bins (core#1171)**: Dart `Conditional color coding: Bin... accepts only 1 to 500 bins`.
- **A value on a shared bin edge, a Min typed above Max, a cycle of linked colourings, a linked source renamed or
  removed (core#1139), stop values typed out of order**: colour arithmetic, not UI checks.
