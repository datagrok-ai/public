---
feature: tileviewer
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas:
  - tileviewer.cp.tiles-font-applied-to-rendering
  - tileviewer.cp.auto-generate-on-columns-change
  - tileviewer.cp.table-rebind-regenerates-form
  - tileviewer.cp.context-menu-inventory
  - tileviewer.cp.viewer-local-filter-vs-dataframe-filter
  - tileviewer.cp.scroll-survives-added-viewer
  - tileviewer.int.viewer-local-filter-vs-df-filter
realizes:
  - viewers.tile-viewer
realized_as:
  - tile-viewer-spec.ts
related_bugs:
  - id: GROK-20096
    status: fixed
  - id: GROK-14215
    status: fixed
  - id: GROK-18230
    status: fixed
  - id: GROK-20207
    status: fixed
  - id: GROK-17775
    status: fixed
  - id: GROK-20376
    status: fixed
  - id: GROK-19950
    status: fixed
  - id: GROK-19983
    status: fixed
  - id: GROK-19016
    status: fixed
  - id: GROK-17081
    status: fixed
  - id: GROK-16598
    status: fixed
---

# Tile Viewer tests

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add Tile Viewer

## Default form rendering

1. Click the first visible tile — row becomes current in the linked grid
2. Click the second visible tile — current row updates to that row

## Row selection

1. Click the first visible tile — tile is highlighted as current
2. Shift-click the third visible tile — tile is highlighted as selected
3. Ctrl-click the fifth visible tile — tile 5 is added to the selection 


## Lanes

1. Open the property panel for Tile Viewer
2. In Property Panel > Data > Lanes, select RACE — four lanes appear (Asian, Black, Caucasian, Other), each labeled by its race value
3. Verify four lanes appear (Asian, Black, Caucasian, Other), each labeled by its race value
4. In Property Panel > Data > Lanes, select SEX — grouping updates to two lanes
5. In Property Panel > Data > Lanes, clear the selection — all tiles return to a single flat lane

## Row source

1. Open the property panel for Tile Viewer
2. In Property Panel > Data > Row Source, set Row Source to Selected
3. Select rows 1–5 in the linked grid — tiles show only those rows
4. In Property Panel > Data > Row Source, set Row Source to Filtered
5. Apply a filter: SEX = M — tiles update to show only matching rows
6. In Property Panel > Data > Row Source, set Row Source to All — all tiles are shown
7. Remove the filter

## Tiles font

1. Open the property panel for Tile Viewer
2. In Property Panel > Style > Tiles Font, change the font size to 18px — tile text and lane headers grow
3. In Property Panel > Style > Tiles Font, change the font family to Arial — the lane headers and the tile text both update to the new font
4. In Property Panel > Style > Tiles Font, reset to the default value (normal normal 13px "Roboto")

## Auto-generate on column change

1. Note which columns the card shows: a freshly added Tile Viewer builds the card by itself (Auto Generate is on), and the card fits at most 10 columns, so on demog (11 columns) one column has no field on the tiles
2. Remove one of the columns that IS on the card: right-click that field on a tile and choose Remove
3. Verify the removed column's field is gone from every tile, and that the column the card had left out has taken the freed slot
4. Open the viewer's context menu on a tile (focus the tile and press the context-menu key) and choose "Edit Form..." — the form editor opens
5. Select one field on the form canvas, delete it, and click CLOSE AND APPLY — the tiles show the edited card without that field, and Auto Generate is now off
6. Remove another column that still has a field on the card the same way — right-click that field on a tile and choose Remove — the field of the removed column goes away here too (a field with no column behind it cannot be drawn), but nothing takes the freed slot and the remaining fields keep their composition and order: the designed layout is not rebuilt

## Multiple table switching

1. Open a second dataset (System:AppData/Chem/tests/spgi-100.csv, table name 'spgi-100')
2. Open the property panel for Tile Viewer
3. In Property Panel > Data > Table, select spgi-100 — tiles update to display spgi-100 rows
4. In Property Panel > Data > Table, select demog — tiles show demog rows again

Switching the table always rebuilds the card for the new table's columns, whether Auto Generate is on or off, so this section checks the rebind and is not a check of auto-generation.

## Context menu items

1. Open the viewer's context menu on a tile: focus the tile and press the context-menu key (a right-click opens the field menu instead)
2. Verify "Edit Form..." item is present (one of the viewer's three own context-menu contributions)
3. Verify "Lanes" item is present (the lanes-column property item, in its own group)
4. Verify "Show Empty Lanes" item is present

## Viewer title and description

1. Open the property panel for Tile Viewer
2. In Property Panel > Description > Title, enter "Patient Cards" — title appears in the viewer header
3. In Property Panel > Description > Description, enter "Demographic data per patient" — description appears below the title
4. In Property Panel > Description > Description Position, select Bottom — description moves to the bottom of the viewer
5. In Property Panel > Description > Title and Property Panel > Description > Description, clear both values — the header is left without a title and the description text disappears

## Filter interaction

1. Open the filter panel (Ctrl+F)
2. Add a filter: SEX = M
3. Verify tiles update to show only male patients
4. Add a second filter: AGE > 50
5. Verify tiles are further reduced
6. Remove all filters — all tiles are restored
7. Close the filter panel

## Viewer filter formula

1. Note the ages shown on the cards — some patients are 50 or younger
2. Open the property panel for Tile Viewer
3. In Property Panel > Data > Filter, enter the formula `${AGE} > 50` and commit it — every card now shows a patient older than 50
4. Verify the grid is unaffected: it still shows every row, because this filter belongs to the viewer and not to the table
5. In Property Panel > Data > Filter, clear the formula — the younger patients are back on the cards

## Scroll position when another viewer is added

1. Make sure the tiles are in a single lane (clear the Lanes grouping)
2. Scroll the lane down with the mouse wheel until it is well away from the top, and note which patient is the first card shown
3. Add another viewer to the same view (for example a histogram)
4. Verify the lane is still scrolled to the same place and still starts with the same patient — the tiles are not re-based to the top

## Automation notes

- PRODUCT FACTS ARE NOT REPEATED HERE. Which categories ship collapsed, why a row query must
  not be scoped inside its category, the font families on offer and what the header does with
  an empty title live in the refdoc `.claude/skills/grok-browser/references/viewers/tile.md`
  ("Property grid (F4 / Context Panel)" and the title-bar table under "Properties"); this
  section carries test decisions only and cites refdoc SECTION NAMES.
- DRIVEN FOR REAL — every "Open the property panel for Tile Viewer" step presses the viewer's
  gear (`clickViewerTitlebarIcon`) and expands the named category (`ensurePropertyCategory`):
  a collapsed category's rows are unreachable. Press the gear directly, never
  `openViewerProperties`, which returns early whenever any property grid is on screen —
  including a stale one, which is exactly the case that needs a reopen.
- DRIVEN FOR REAL — Row Source, Table, Title, Description, Description Position and Tiles Font
  through their own property-grid editors, so a broken or uncommitted editor fails the step
  instead of passing silently. Tiles Font goes through the font widget's size input and family
  select; the family step picks Arial because that select's whole offer is the four families
  the refdoc lists for `prop-tiles-font`.
- PROCEDURE — reach a Data row by waiting for the row itself, not through the expander; expand
  collapsed categories first; drive a row only once it is really VISIBLE (press the gear again,
  re-expand the category).
- PROCEDURE — a step that acts on the table does so BEFORE it reaches for a property row:
  acting on the table rebuilds the Context Panel (refdoc: "Acting on a property REBUILDS the
  Context Panel").
- PROCEDURE — both panels are shared state, so a step brings them to the state it needs instead
  of assuming it. Create the filter panel only when it is absent, and never wait for a filter
  widget before adding one: in the later sections the panel is already on screen with an empty
  filter group, so that wait never returns.
- The honest expectation after clearing a title through the panel is a header with NO title at
  all — not a return to the type name.
- SUBSTITUTED ACTUATION, with the reason each surface is unavailable:
  - Lanes (Property Panel > Data > Lanes) — not addressable from the DOM (refdoc: the
    `prop-lanes` row under "Property grid (F4 / Context Panel)"), so the property is set
    through the viewer's JS API (`viewer.props.lanesColumnName`); Lanes steps 2, 4 and 5 are
    'RACE', 'SEX' and null. ASSERT = the lane structure the product rendered, never the
    property value.
  - the individual filters inside the filter panel — added through the filter group, because
    the categorical and histogram widgets are canvas controls with no per-category selector.
    The panel itself is real: `openFilterPanel` stands in for Ctrl+F and waits for a real
    filter to render, and it is closed through the panel's own title-bar close button.
  - selecting rows 1–5 for the Row source section — set on the DataFrame's selection, because
    the grid is canvas with no per-row selector. ASSERT = the tile population the product
    rendered.
- COLUMN-REMOVING GESTURE IS REAL: right-click the field's value on a tile and pick the removal
  leaf the refdoc names under "Field (column) menu — what a real right-click opens". A
  programmatic removal reaches the same downstream listener while leaving the control
  unexercised.
- KNOWN GAPS: the context-menu path is covered only for the viewer's own three contributions in
  the viewer menu — the field menu's own inventory and the rest of the viewer menu are not
  enumerated here. The two filter-axis steps assert the tiles' CONTENTS, not their number: the
  lane is a virtualized window, so a tile count says nothing about how many rows passed the
  filter.

Spec must keep:

* An auto-generated card holds at most 10 fields, picked by relevance — on demog (11 columns)
  one column (SEVERITY) gets no field, and the card shows 10 labels plus 10 values.
* Never use "add a column and expect it on the card" as the auto-generate signal: the cap is
  already used up, so the form regenerates (field order changes) without selecting the new one.
* The usable auto-generate signal is deleting a column that IS on the card: its field leaves
  every tile and the freed slot refills from a column that had been left out.
* Read which column is on the card from the card itself, never assume it by name — otherwise the step addresses a column the card does not show.
* Pick the delete victim from the rendered fields but never simply the first: a boolean column
  renders as a disabled checkbox whose field-menu click times out. Filter the CHOICE of victim
  only — the composition read must still count every field.
* Scope a card's field reads to the root of the viewer under test: the page-wide viewer selector
  resolves to the first Tile Viewer on the page (the main one on demog), so a check meant for a
  second viewer finds zero fields there and passes or fails for the wrong reason.
* What tells the two states apart is the REFILL, not the survival of the deleted column's field:
  a field whose column is gone disappears in both states.
* Never assert that the designed state's field set is UNCHANGED — too strong, it will fail. The
  designed-state bar is that no left-out column takes the freed slot and the survivors keep
  their composition and order.
* Right-clicking a tile opens the FIELD menu (no Edit Form entry), so a wait for the viewer menu
  after a right-click never returns.
* The viewer menu opens by focusing a tile's sketch and pressing the context-menu key, or from
  the viewer hamburger — which lives in the panel next to the viewer root, not inside it.

---
{
  "order": 18,
  "datasets": ["System:DemoFiles/demog.csv"]
}
