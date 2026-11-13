# Box plot — what the features cannot say yet

Written during the bar chart / box plot / histogram / pie / line chart / density plot gap round of
23 Oct 2026, against the box plot's "Automation surface" (`core/client/d4/lib/src/viewers/box_plot/CLAUDE.md`)
and `files/TestTrack/Viewers/BoxPlot/box-plot-ui.md`; trimmed in the BDD PR round of 8 Nov 2026.

## Resolved

- The violin style over a value column holding a raw NaN (core#960): ApiTests `AI: Viewers: BoxPlot extras`.

## Still open

- **The mouse-over highlight follows the pointer** (`box-plot-ui.md` "Mouse-over highlight"): needs a second marker area
  (`marker of row <r>`) or a `hovered row` reading from the box plot's mouse-over hit test in `box_plot_core.dart`.
- **The other groups fade under a cross-viewer hover** (github-3065): the fade is a translucent overlay the colour step
  cannot tell apart; needs a `faded groups` reading in the box plot's status.
- **Each category's statistics under its own box, in the column's format** (public#3693, public#3773): needs
  `stats <category>` areas and a reading of the text drawn in each.
- **The log-axis warning** ("… can't be visualized in the log scale") has no area; needs a `y warning` area, as the
  density plot reports.
- **The info mark of a one-row group under Show Group Comparison** (decided 2026-10-28; its tooltip says the group has
  only 1 value): needs `marker of <group>` areas (`_markerHits` in `box_plot_group_comparison.dart`).

## Won't do

- A plot shrunk under Auto Layout and hovered (core#960): ApiTests `AI: Viewers: BoxPlot extras` "a pointer over a box
  plot shrunk under Auto Layout logs no error (core#960)".
- Alpha typed in the panel, Category 1 removed, every row filtered out, listener counts: the shared viewer rules (Dart
  `Data Manipulation`, `Lifecycle`). The control group filtered out, the comparison at 30 and 31 categories, the
  comparison under a filter: observations, not gaps.
