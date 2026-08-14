---
feature: grid
target_layer: playwright
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: [viewers.grid, grid.sort, grid.column-resizing, grid.column-reordering, grid.row-selection, grid.column-selection, grid.editors, grid.copy-paste, grid.popup-menu, grid.column-tools, grid.keyboard-navigation, grid.special-rows, powergrid.action.sparklines, powergrid.cell.sparklines, viewers.scatter-plot]
realized_as:
- grid-spec.ts
related_bugs:
- id: GROK-18256
  status: fixed
- id: GROK-19942
  status: fixed
- id: GROK-19717
  status: fixed
expected_results:
- anchor: Column Sizing
  expectation: Choosing Optimal from the Column Sizing submenu resizes every column to fit its content,
    and Minimal and Maximal each change the widths again in their own direction.
- anchor: Header Histogram
  expectation: Adding Top > Histogram puts a miniature histogram strip in the column header area; repeating
    the path removes it.
- anchor: GROK-18256
  expectation: Removing a summary column with the top-panel icon takes it out of the grid — the GROK-18256
    guard, where only the menu path worked.
- anchor: GROK-19942
  expectation: After the column a Tags summary cell refers to is removed, the grid keeps rendering and
    raises no error — the GROK-19942 guard, where the whole grid stopped drawing.
- anchor: GROK-19717
  expectation: A project holding a table of extracted rows reopens with its grid and no error balloon
    — the GROK-19717 guard.
- anchor: Added Viewer
  expectation: A second Grid added from the toolbox shows the table's row count, takes its own row-height
    and column-label settings without touching the view's own grid, re-binds to spgi-100 when Data >
    Table is switched, and leaves the view's grid intact when closed.
- anchor: Multi-Column + Row-Height Resize
  expectation: With two columns selected, dragging one border resizes both; dragging the boundary between
    two row headers changes the row height; dragging a border far left shrinks the column until its
    values stop being legible.
- anchor: Tooltip Settings
  expectation: The Tooltip submenu switches the column tooltip between Default, None, Form and Columns,
    and the Columns choice stores the chosen column list and shows exactly those values on hover.
---

# Grid tests (Playwright)

## Setup

1. Close all open views.
2. Open demog.csv. Every section below starts from this state unless it names a different table,
   in which case its own first step opens it.


## Column Sizing from the Context Menu

1. Right-click a data cell, expand Column Sizing, and choose Optimal — every column resizes to fit
   its content
2. Repeat with Minimal, then Maximal — the widths follow each choice

## Header Histogram Strip

1. Right-click a data cell, expand Add > Top, and choose Histogram — a miniature histogram strip
   appears in the column header area
2. Repeat the path and turn it off — the strip disappears

## Removing a Summary Column by the Top-Panel Icon

*Regression guard for GROK-18256, where this removal path did nothing while the menu path worked.*

1. Right-click a data cell, expand Add > Summary Columns, and add Sparklines
2. Select the new summary column and remove it with the remove icon on the top panel — the column
   disappears from the grid

## Summary Column Surviving a Removed Source Column

*Regression guard for GROK-19942, where the whole grid stopped rendering.*

1. Right-click a data cell, expand Add > Summary Columns, and add Tags
2. Right-click the CONTROL column header and choose Remove — the Tags cell was referring to it
3. Scroll the grid — it keeps rendering and raises no error

## Extracted Rows in a Saved Project

*Regression guard for GROK-19717, where such a project failed to reopen.*

1. Ctrl+click the row-number cells of rows 1 through 5 so each joins the selection
2. Use Select > Extract Selected Rows — a new table opens
3. Save the view as a project with the ribbon Save button, then Close All
4. Reopen the project — it opens with its grid and no error balloon

## Grid as an Added Viewer

*The rest of this section's coverage acts on the table view's own grid. This subsection is the
only place a SECOND, explicitly added Grid viewer is exercised — it has a narrower surface (its
column menu is off by default) and is not a substitute for the main grid.*

1. Add a Grid from the toolbox Viewers section — a second grid appears alongside the first and
   shows the same number of rows as the table
2. Open the added viewer's settings with its own gear icon and change Row Height — only the added
   viewer's rows change height
3. Toggle Show Column Labels off and back on in the same panel — the added viewer's headers hide
   and reappear
4. Open spgi-100.csv so a second table is available, go back to the demog view, then open the
   Context Panel for the added viewer and switch Data > Table to spgi-100 — the viewer re-binds
   and shows spgi-100's row count while the view's own grid still shows demog
5. Close the added viewer — the table view keeps its own grid

## Multi-Column and Row-Height Resizing by Drag

*Dataset: demog*

1. Put AGE and HEIGHT into the column selection — selecting columns BY MOUSE is a manual case
   (see the manual checklist), so a test may establish the selection directly — then drag the right
   border of the AGE header 40 pixels to the
   right — both selected columns end up at the same new width, and the columns around them keep
   the width they had
2. Press Esc to drop the column selection, then drag the AGE border again — this time only AGE
   changes width
3. Drag the border between the first two row numbers 30 pixels downwards — every row grows taller
   by that amount; the column widths stay exactly as they were
4. Drag the AGE right border far to the left until the column is only a few pixels wide — the
   values stop being drawn as numbers and appear as small circles instead
5. Drag the WEIGHT right border left as far as the column's own left edge — the column collapses
   to a hairline showing no content, while still counting as a visible column

## Column Tooltip Settings

*Dataset: demog*

1. Right-click the AGE column header and open Tooltip > Current Column — the choices are Default,
   Form, Columns and None, and the one in force is marked
2. Choose None, then hover a cell in the AGE column — no tooltip is shown
3. Right-click the AGE header again, open Tooltip > Current Column and choose Columns — a
   "Select columns..." dialog opens
4. Press All in that dialog and confirm with OK — the dialog closes
5. Hover a cell in the AGE column — the tooltip now lists every chosen column with this row's
   values
6. Right-click the AGE header, open Tooltip > Current Column and choose Default — the mark moves
   back to Default



## Automation notes

- The hamburger inline filter went back to grid-ui.md. Recon drives it — the handles are painted on
  the histogram's bottom strip and a drag there does move the filtered row count — but a spec could
  not reproduce that grab reliably. Two of the three sections recon recovered stayed here; this one
  did not.

- The spec must cover THIS file's sections and nothing else. Steps left over from the sections that
  moved into focused scenarios — row and column selection, cell editing, clipboard, both context
  menus, cell renderers, frozen-column properties, the hamburger stats popup and search — are to be
  deleted from the spec, not kept "just in case": each of them is already driven by a focused
  scenario, and a second copy here means two places to fix when the UI moves.

- Three sections returned here from the manual checklist after a live recon refuted their old
  "not automatable" reasons: linked multi-column resize with the row-height drag, the tooltip
  settings, and the hamburger inline filter. Two others were refuted as behaviour but stayed manual
  because a driver could not reach them — the form designer and the column permissions — and they
  live in grid-ui.md, not here.
- Clearing a column selection with the property alone leaves the grid's resize LINK in place: a
  later border drag still resizes the previously selected columns. Only a trusted Escape clears it.
- The sections this file used to carry for row/column selection, cell editing, clipboard, keyboard
  navigation, the two context menus, cell renderers, frozen-column properties, the hamburger menu
  and search were removed once their focused scenarios landed — this file must not restate what a
  focused scenario already drives. What stays here is the residue no focused scenario owns.


- Two of the sections here came back from the manual checklist after a live recon refuted their old
  "not automatable" reasons: the linked multi-column resize with the row-height drag, and the
  tooltip settings. The form designer, the column permissions and the hamburger inline filter were
  refuted as BEHAVIOUR too, but stayed manual because a driver could not reach them — they live in
  grid-ui.md and their sections there say exactly which path failed.


- Cross-viewer row-state synchronisation was dropped deliberately: current row and selection are
  dataframe state shared by every viewer, so there is nothing grid-specific to prove.

- grid-spec.ts still carries the legacy end-to-end walkthrough written against the pre-split
  sections and must be regenerated for the shape above. Its import depth was repaired
  ('../../spec-login', '../../helpers/viewers') — the nested-section form — so it can at least be
  collected and run.

---
{
  "order": 15,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
