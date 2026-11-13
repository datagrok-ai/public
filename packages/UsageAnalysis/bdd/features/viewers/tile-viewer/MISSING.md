# Tile viewer — what the features cannot say yet

Written during the viewer-gaps round of 23 Oct 2026, against
`files/TestTrack/Viewers/TileViewer/tile-viewer-ui.md` and the viewer's "Automation surface"
(`core/client/d4/lib/src/viewers/tile_viewer/CLAUDE.md`). New this round:
`tile-viewer-designed-rebind.feature` (a card designed on spgi-100, then Table → demog-1000) and
`tile-viewer-lane-scroll.feature` (GROK-20912, `@known-failure`).

## Still open

- Colour coding on tiles: a card field takes its column's colour as the input's background; needs a `background of
  <COL> of row <r>` reading, as the forms viewer's `background of <COL> of <label>`.
- core#1151 — a value the card cuts off shows no tooltip (public#1593, the value half; the caption half is claimed in
  `tile-viewer-mirroring.feature`).
- A card dragged between lanes while a tile viewer is open in another view logs `NullError: method not found: 'cZ' on
  null` (2026-11-03): to be filed.

## Won't do

- Lanes picked in the Context Panel: the property grid's column picker is generic, claimed on the scatter plot
  (`sunburst.feature`) and the trellis plot (`trellis-plot-viewer-filter-and-menus.feature`).
- The lane of an empty value titled "(no value)" (core#1217): xamgle `d4 | Viewers | Fixed issues | Tile viewer titles the
  lane of an empty value`.
- Escape in the form designer: a view, not a dialog; Escape changing nothing is no claim.
- A removed or renamed lanes column, Row Source Selected with nothing selected, a one-row table: the rules every viewer
  shares (`../MISSING.md`).
