# Trellis plot — what is still missing

The gap round of 2026-10-07 over TestTrack `Viewers/TrellisPlot/` (six md cases, their specs and
`trellis-plot-ui.md`) against the nine features already in this folder. Its features are the seven
added that day: `use-in-trellis`, `viewer-filter-and-menus`, `selectors-and-full-screen`,
`curves-table`, `click-gaps`, `scales-and-paging`, `pick-up-color-and-scroll`. Triaged with the lead on
2026-10-11: what is not planned is under "Won't do", with the reason.

## Still open

- Show Empty Bars off gives every inner bar chart its own categories (public#3616): the trellis reports a cell's picture,
  not its bars; needs a per-cell `bars of <cell>` reading from the inner bar chart's status.

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
- **core#1057's edge cases** (a label click under a filter, the arrow keys after a table switch, Tiles Per Row 0, an
  inner viewer type in another case, the leaks): xamgle `d4 | Viewers | Trellis plot | …` and `d4 | Viewers | Table swap |
  Trellis plot…`; Tiles Per Row past its bounds: core#1200's clamp tests (`../MISSING.md`).

## Traps found

- A cell already current is not clicked again by a scenario that expects a selection: pick another cell.
- Pick Up / Apply, inner color, scroll: `x label <category>` / `y label <category>` areas are the way to
  tell which categories a scroll brought into the window.
