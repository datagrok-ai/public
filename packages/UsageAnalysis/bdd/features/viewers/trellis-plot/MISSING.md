# Trellis plot — what is still missing

The gap round of 2026-10-07 over TestTrack `Viewers/TrellisPlot/` (six md cases, their specs and
`trellis-plot-ui.md`) against the nine features already in this folder. Its features are the seven
added that day: `use-in-trellis`, `viewer-filter-and-menus`, `selectors-and-full-screen`,
`curves-table`, `click-gaps`, `scales-and-paging`, `pick-up-color-and-scroll`. Triaged with the lead on
2026-10-11: what is not planned is under "Won't do", with the reason.

## Resolved

- **Reset Inner Range Sliders is Reset View** since core GROK-17986 (2026-10-08): the trellis's menu has one Reset
  View that sets the category scrolls and every inner range slider back; `scales-and-paging` picks it (review of
  2026-10-11).
- **The table and the viewer's Filter formula are set in the context panel**, as a person does, in `curves-table`
  and `viewer-filter-and-menus` (review of 2026-10-11).
- **Cell clicks after On Click and Row Source were picked in the context panel** — not reproduced on 2026-10-11
  (core master `6c4378c593`): after `click-gaps`'s context-panel sequence, On Click = Select and a click on
  M | Asian select its 8 rows.
- **Escape after the inner viewer changed — not a defect** (walked by hand on dev, 2026-10-07). Changed in the
  viewer selector, as a person does, the current cell and the filter stay and Escape drops the filter. Only a
  write of the Viewer Type property (API, or a step that sets the property) resets the current cell, after
  which the Escape handler (`trellis_plot_core.dart`, `if (currentCellPos != null)`) has nothing to act on;
  `click-gaps` changes the type through the selector.
- **Whether a selector strip is visible** (2026-10-12): the `x selectors` / `y selectors` and `x/y selector <n>`
  areas are reported only while the strip is on the page (`trellis_status.dart`, `htmlGetVisible`), not from the
  layout's flag. `selectors-and-full-screen` claims Selectors step 5: Show X Selectors off, Auto Layout's shrink
  and restore leave the X strip off, the property brings it back.
- **What the full-screen dialog shows** (2026-10-12): the `full screen rows` reading is the row count of the inner
  viewer the dialog shows, present only while it is open; `selectors-and-full-screen` claims 480, then its absence
  after Escape.
- **The rows of the table a trellis is bound to** (2026-10-12): the `table rows` reading; `curves-table` claims 30
  after the switch to curves and 1000 after the switch back.

## Still open

Nothing planned.

## Won't do

- **The (+) column picker's hover preview, Escape and blank entry** (GROK-19673, `trellis-plot-split-and-pick-inner.md`
  section 1 steps 5-7): needs a picker-opened signal on the (+) and a hover step inside an open column picker — a
  core signal and a library step for one case; the split itself is claimed through the selectors.
- **Pick Up / Apply step 8** (a range slider moved on the second trellis leaving the first alone): needs the area
  steps to resolve hover-revealed sliders of a second trellis; the copy itself is claimed.
- **The Row Source ladder through the property list** (`trellis-plot-click-to-filter.md` section 2 step 5,
  GROK-13205): `trellis-plot-row-source` walks all eight rungs through the API; `click-gaps` sets Row Source and
  On Click in the panel for the correction between them.
- **A (+) click at the end of its range**: a "nothing happened" claim has no signal to wait for; `scales-and-paging`
  claims the icon's `aria-disabled` instead.
- **A floating viewer after a layout is applied, undocking, browser zoom** (`trellis-plot-ui.md`): no step undocks
  a viewer into a floating window or zooms the page, for one checklist line.
- **The ribbon Save and the Layout menu** (`trellis-plot.md` "Layout and Project save/restore"): the persistence
  features save through the API, and what a layout stores is an API test's claim; applying a saved layout from the
  Layout menu is the layouts features' subject.
- **Pie chart "Marker Color"** (`trellis-plot-ui.md` "Inner viewer color coding"): the pie chart has no such
  property; its slices are coloured by its Category, which `pick-up-color-and-scroll` sets.
- **Multi Curve steps 5-7** (curve X/Y, paging, zoom slider inside the cells): the case itself leaves them manual
  for want of a recon of the curve viewer's controls.

## Traps found

- A cell already current is not clicked again by a scenario that expects a selection: pick another cell.
- Pick Up / Apply, inner color, scroll: `x label <category>` / `y label <category>` areas are the way to
  tell which categories a scroll brought into the window.
