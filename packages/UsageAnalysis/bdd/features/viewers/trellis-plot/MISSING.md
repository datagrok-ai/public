# Trellis plot — what is still missing

The gap round of 2026-10-07 over TestTrack `Viewers/TrellisPlot/` (six md cases, their specs and
`trellis-plot-ui.md`) against the nine features already in this folder. Its features are the seven
added that day: `use-in-trellis`, `viewer-filter-and-menus`, `selectors-and-full-screen`,
`curves-table`, `click-gaps`, `scales-and-paging`, `pick-up-color-and-scroll`. Written with today's
steps only; what could not be said is listed here.

## Resolved

- **Escape after the inner viewer changed — not a defect** (walked by hand on dev, 2026-10-07). Changed in the
  viewer selector, as a person does, the current cell and the filter stay and Escape drops the filter. Only a
  write of the Viewer Type property (API, or a step that sets the property) resets the current cell, after
  which the Escape handler (`trellis_plot_core.dart`, `if (currentCellPos != null)`) has nothing to act on;
  `click-gaps` changes the type through the selector.

## Missing steps and signals

- **The (+) column picker** (GROK-19673, `trellis-plot-split-and-pick-inner.md` section 1 steps 5-7): its
  hover preview, Escape undoing the preview, the click that commits, and the blank entry that drops a
  split column. The picker does not open reliably from a pointer gesture on the `add x column` area
  (see `trellis-plot-split-and-inner-type`), and no step hovers a row of an open picker and then leaves
  with Escape. Wanted: a picker-opened signal on the trellis's (+) and a "hovers over {string} in the open
  column picker" step.
- **Pick Up / Apply step 8**: a range slider moved on the second trellis leaving the first alone. It needs
  Global Scale and a hovered cell on the second trellis while the first is in the view; the area
  steps address the first matching viewer's slider. Wanted: the area steps resolving `second trellis
  plot viewer` for hover-revealed sliders (or a reading of each trellis's slider range).
- **Whether a selector strip is visible** (`trellis-plot.md` "Selectors", and step 5, the X strip that
  stays off through Auto Layout's shrink and restore, which reads only the flag and so is not claimed): the `x selectors` / `y selectors`
  areas follow the layout's own `_showXSelectors` / `_showYSelectors` flag (`trellis_status.dart`), not
  whether the strip is visible on the page, and the `x selector <n>` areas are reported even while hidden.
  Wanted: both gated on `htmlGetVisible`, or an `x selectors shown` reading of the pickers a person can see.
- **What the full-screen dialog shows**: the dialog's title names the cell (`SEX: F, RACE: Caucasian`) and its
  canvas is painted, but the viewer inside reports no rows; wanted: a reading of the dialog viewer's rows.
- **The rows of the table a trellis is bound to**: `should be bound to table` reads the table's name; the case
  also asks for its row count after the switch to curves. Wanted: a `table rows` reading.
- **The Row Source ladder through the property list** (`trellis-plot-click-to-filter.md` section 2 step 5,
  GROK-13205): `trellis-plot-row-source` walks all eight rungs through the API; `click-gaps` sets Row Source
  and On Click in the panel only for the correction between them.
- **A (+) click at the end of its range**: the click has nothing to add and no step waits for a "nothing
  happened" signal, so `scales-and-paging` claims the icon's `aria-disabled` state and makes no such click.
- **A floating viewer after a layout is applied, undocking, browser zoom** (`trellis-plot-ui.md`): no step
  undocks a viewer into a floating window or zooms the page.
- **The ribbon Save and the Layout menu** (`trellis-plot.md` "Layout and Project save/restore"): the
  persistence features save through the API. The project Save dialog has steps; applying a saved layout
  from the Layout menu has none.

## Not translated, by choice

- **Pie chart "Marker Color"** (`trellis-plot-ui.md` "Inner viewer color coding"): the pie chart has no such
  property; its slices are coloured by its Category, which `pick-up-color-and-scroll` sets.
- **Multi Curve steps 5-7** (curve X/Y, paging, zoom slider inside the cells): the case itself leaves them
  manual for want of a recon of the curve viewer's controls.

## To check by hand

- **Cell clicks after On Click and Row Source were picked in the context panel.** After the context-panel
  scenario of `click-gaps` (Row Source = Filtered, On Click = Filter, a cell click, Escape, Row Source =
  Filtered again, which moves On Click to None), a later On Click = Select and a click on M | Asian left no
  current cell and no selection. The same sequence through the API leaves later clicks working. Whether a
  person sees it too is unknown; `click-gaps` puts its context-panel scenario last.

## Traps found this round

- A cell already current is not clicked again by a scenario that expects a selection: pick another cell.
- Pick Up / Apply, inner color, scroll: `x label <category>` / `y label <category>` areas are the way to
  tell which categories a scroll brought into the window.
