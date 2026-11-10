# Tooltips — what is still missing

The features translated from TestTrack `Tooltips/` (seven md cases and their specs, six features:
uniform-default-tooltip is claimed by default-tooltip-visibility, whose first scenario hovers the same three
viewers with the same setup). Run on the local master stand, 2026-10-07.

## Resolved

- The grid's own default tooltip shows nothing — by design (2026-10-07): the features switch the grid to "inherit from
  table" first, and grid-visible-columns-tooltip claims the default.
- The table's tooltip lists at most seven columns — by design (2026-10-27).
- Tooltip > Hide switches off every tooltip of the table — by design (2026-10-07).
- A custom tooltip with an empty Row Tooltip lists only the viewer's own data columns — by design (2026-10-07).
- The scatter plot lists its axis columns first under Data Values = Merge — by design (2026-10-11).

## Won't do

- The trellis plot's own tooltip after Hide: its inner viewers' marks are not areas of the trellis; the line chart shows
  no tooltip until Edit Aggregated Tooltip configures one.
- The Edit Tooltip dialog's list has three columns (type, name, box): a reading of a dialog's inner grid for one low-value
  line.
- default-tooltip step 4, widening the context panel: grid-visible-columns-tooltip pushes the columns off by the column
  resizer.
- line-chart-aggregated-tooltip's "in-viewer filter": named in the case's title only.
- What the tooltip's renderer draws (molecules in SPGI's tooltip): a renderer, not a feature.
- The `.tooltip` and `.tooltip-visibility` tags the old specs read: echoes of what the gesture wrote.
- CANCEL in Edit Tooltip, every column unticked, the Add aggregation default: low value (a dialog writes on OK only).
- A renamed or removed tooltip column (core#1215): xamgle `d4 | Tooltip | Tooltip columns follow a column rename and skip
  a removed one`.
