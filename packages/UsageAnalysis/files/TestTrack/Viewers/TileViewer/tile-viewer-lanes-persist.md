---
feature: tileviewer
realizes_atlas:
  - tileviewer.cp.lanes-cells-and-layout-persist
  - tileviewer.int.lane-drag-writes-dataframe-cell
  - tileviewer.int.lanes-rebuild-vs-restyle-scope
  - tileviewer.int.column-rename-rewrites-sketch-state
realizes:
  - viewers.tile-viewer
priority: p0
target_layer: playwright
coverage_type: smoke
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
realized_as:
  - tile-viewer-lanes-persist-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 2"
    expectation: >-
      Exactly one lane exists carrying class .d4-tile-viewer-lane-single; no
      .d4-tile-viewer-lane-header is visible (hidden by CSS).
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      Four lanes appear, each classed .d4-tile-viewer-lane-multi; the ordered
      .d4-tile-viewer-lane-header texts are exactly ['Asian', 'Black',
      'Caucasian', 'Other'] matching the RACE category order; no
      .d4-tile-viewer-lane-single class is present.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      The ordered lane header texts equal the explicit list ['Black', 'Asian'] —
      NOT the column's category order — and only those two lanes exist;
      'Caucasian' and 'Other' have no lane.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      With RACE filtered to a single category (Black) and showEmptyLanes=true
      (default), both explicitly listed lanes ['Black', 'Asian'] remain present;
      the filtered-out 'Asian' lane holds zero .d4-tile-viewer-form tiles while
      the unfiltered 'Black' lane holds its tiles.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      After setting showEmptyLanes=false, only the lane(s) with rows survive;
      the empty lanes are absent from the DOM. After restoring
      showEmptyLanes=true the empty lanes return in their original order
      (GROK-20096 regression guard).
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      After clearing lanesColumnName, exactly one lane exists carrying
      .d4-tile-viewer-lane-single; the filtered row set is shown in that single
      lane.
  - anchor: "Scenario 2 Step 2"
    expectation: >-
      After dragging a tile from lane A to lane B, the dragged row's RACE cell
      in the DataFrame equals lane B's category value, and df.currentRowIdx is
      that row; the tile now appears in lane B's .d4-tile-viewer-lane-content
      and is absent from lane A's.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Dragging a tile onto its own current lane leaves the RACE cell unchanged —
      the cell value equals the lane it was already in.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      With allowDragBetweenLanes=false, performing the same drag gesture leaves
      the RACE cell unchanged.
  - anchor: "Scenario 3 Step 2"
    expectation: >-
      After editing the AGE cell in the grid, the tile's
      input[name="input-AGE"].value equals the grid's displayed cell text for
      the same row (the formatted string, not the raw numeric value) without
      reopening Edit Form (GROK-17775 regression guard).
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      After renaming AGE to AGE_YRS, the tile's label input
      (input.d4-sketch-column-name) reads 'AGE_YRS', the value input is now at
      [name="input-AGE-YRS"] and still resolves to the same column's display
      string, and the host attribute is [name="div-AGE-YRS"] (GROK-20207
      regression guard). A tile for a column NOT renamed is unchanged.
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      After the tile form is brought to its 10-field cap so the calculated float
      column is promoted onto the card, its tile input[name="input-COMPUTED-H"]
      (underscore slugged to hyphen) .value equals the grid's displayed cell text
      for the same row and differs from the raw cell value (GROK-20376 regression
      guard); the set of fields shown on the tile has changed from the baseline.
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      After re-applying the saved layout, the Tile Viewer is present with the
      lanes column set and the lane structure re-rendered (ordered
      .d4-tile-viewer-lane-header texts match the saved configuration); the
      console is clean — no 'method not found' or null-method errors (GROK-18230
      / GROK-18037 regression guard).
  - anchor: "Scenario 4 Step 6"
    expectation: >-
      After project reopen, the Tile Viewer shows the same lane structure and a
      spot-checked tile's value input equals the column's display string for
      that row.
---

# Tile Viewer — Lanes ladder, tile-content mirroring, and persistence

## Setup

1. Close All open views and tables.
2. Open System:DemoFiles/demog.csv (11 columns: USUBJID, AGE, SEX, RACE, DIS_POP,
   HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY; RACE has 4 categories —
   Asian, Black, Caucasian, Other; SEX has 2; CONTROL is boolean). The
   auto-generated tile form is capped at 10 fields, so one column (SEVERITY) is
   omitted from the card.
3. Add a Tile Viewer to the current table view.
4. Confirm the viewer is in the auto-generated state: no form has been manually
   designed yet (the viewer shows auto-generated tiles for all columns, which is
   the entry condition for Scenario 3 and Scenario 4).

## Scenarios

### Scenario 1: Lanes ladder — structure accumulates through six steps

Steps:
1. Confirm the baseline: no lanes column is set (freshly added viewer shows a
   single lane with all rows).
2. Verify the single-lane baseline — inspect the lane structure in the viewer.
3. Set the **Lanes Column** to **RACE** using the Tile Viewer property panel
   (Properties > Lanes Column, select RACE). Note: the lanes column selector
   opens a canvas column grid with no per-column DOM selector — this step is
   performed programmatically in automation (see Automation notes). Verify the
   four-lane structure once the viewer re-renders.
4. Configure the explicit lanes list to ['Black', 'Asian'] by ordering the lane
   list in the Tile Viewer property panel (Properties > Lanes). Note: the
   explicit lane list has no interactive UI control — this step is performed
   programmatically in automation (see Automation notes). Verify the two-lane
   structure appears in the listed order ['Black', 'Asian'].
5. Apply a categorical RACE filter that keeps a single category ('Black') using
   the registered filter helper (Filter Panel categorical state). With the
   default Show Empty Lanes setting on, verify the empty-lane behaviour on the
   two listed lanes: the filtered-out 'Asian' lane holds zero tiles while 'Black'
   keeps its tiles.
6. In the Tile Viewer property panel, turn off **Show Empty Lanes** (Data > Show
   Empty Lanes toggle — the row lives under the Data category, see refdoc
   "Property grid (F4 / Context Panel)").
   Verify the empty lanes disappear. Then turn **Show Empty Lanes** back on and
   verify the empty lanes return in their original order (GROK-20096 regression
   guard).
7. Clear the **Lanes Column** setting in the Tile Viewer property panel
   (Properties > Lanes Column, clear the selection) — performed programmatically
   in automation, same reason as step 3 (see Automation notes) — and remove the
   RACE filter. Verify the viewer returns to a single lane.

### Scenario 2: Drag between lanes writes the DataFrame cell

Steps:
1. Set the **Lanes Column** to **RACE** in the Tile Viewer property panel
   (Properties > Lanes Column, select RACE). Confirm four lanes appear with
   tiles in each.
2. Identify a specific tile in the 'Asian' lane (note its row index and RACE
   cell value). Drag that tile onto a tile of the 'Black' lane — the lane
   CONTENT area, because the lane header does not accept a drop — using a real
   mouse drag gesture and release. Verify the cell write and tile relocation.
3. Identify a tile in the 'Black' lane. Drag it back onto the 'Black' lane
   (same lane — own-lane drop). Verify the drag has no effect on the cell value.
4. In the Tile Viewer property panel, turn off **Allow Drag Between Lanes**
   (Properties > Misc > Allow Drag Between Lanes — the row lives under the Misc
   category). Attempt to drag a tile from the 'Black' lane to the 'Caucasian'
   lane. Verify the RACE cell is unchanged.
5. In the Tile Viewer property panel, turn **Allow Drag Between Lanes** back on
   and clear the **Lanes Column** setting.

### Scenario 3: Tile content mirrors the DataFrame

Tiles display formatted column values. Each column field on the tile shows the
same display string that appears in the linked grid — not the raw underlying
value. The tile label shows the column name and updates when the column is
renamed.

Steps:
1. Make sure the viewer is in the auto-generated state (verified in Setup step 4).
   Note the displayed value for the AGE field on the tile for row 0 as the
   baseline.
2. In the linked grid, edit the AGE cell for row 0 to a different integer (e.g. 99).
   Without reopening Edit Form, verify the tile's AGE field updates to match the
   grid's displayed cell text for that row.
3. In the linked grid, rename the AGE column to AGE_YRS (double-click the column
   header or use the column menu Rename). Without reopening Edit Form, verify
   the tile label updates to 'AGE_YRS' and the AGE_YRS field still shows the
   same column's display string. Also verify that a column NOT renamed (e.g.
   HEIGHT) has unchanged tile fields.
4. Add a calculated float column named COMPUTED_H with formula `${HEIGHT} * 1.0`.
   The auto-generated form is capped at 10 fields and demog already fills it, so a
   freshly added column is NOT shown until a slot is freed. Bring the tile form
   down to exactly 10 columns (remove two columns) so COMPUTED_H is promoted onto
   the card, then confirm both halves: (a) the set of fields shown on the tile
   has changed from the baseline, and (b) COMPUTED_H's tile value — read from
   input[name="input-COMPUTED-H"] (the underscore is slugged to a hyphen) —
   equals the grid's displayed cell text for row 0 and differs from the raw cell
   value (GROK-20376 regression guard; a raw-value comparison would fail on
   correct behaviour). Remove the calculated column in the teardown. Note: the two
   columns freed to reach the cap (DEMOG and SEVERITY) are NOT restored — their
   data is gone once removed — so the frame stays at nine columns from here on, and
   Scenario 4 saves/persists that nine-column frame.

### Scenario 4: Peak configuration persists through layout save and project save

Precondition (carried from Scenario 3 Step 4): DEMOG and SEVERITY have been
removed and are not restored, so this scenario runs on the resulting nine-column
frame. Lane structure and persistence are independent of those two columns.

Steps:
1. Configure the peak: set the **Lanes Column** to **RACE** in the Tile Viewer
   property panel, then set the explicit lanes list to ['Black', 'Asian',
   'Caucasian'] (see Automation notes for how to set the explicit list), and
   confirm **Show Empty Lanes** is on. Confirm the three-lane structure appears.
2. Save the current view layout through the **View | Layout | Save to Gallery**
   top-menu leaf (driven with the real mouse). Locate the saved layout by diffing
   the table's applicable-layout ids before and after the save.
3. Modify the current view: close the Tile Viewer and add a Grid viewer. The view
   now differs from the saved layout.
4. Re-apply the saved layout (Home > Layouts, select the layout, Apply). Wait for
   the viewer to finish rendering. Verify the lane structure is restored.
5. Delete the probe layout in the teardown (even on failure).
6. From the layout-restored state (Tile Viewer with RACE lanes), save the project
   using the ribbon **Save** button. Note the project name. Close All. Reopen the
   project from the workspace (Browse > Recent Projects or the file browser).
   Spot-check a tile's field value against the grid's displayed cell text for
   that row. Delete the probe project in the teardown.

## Step-to-spec mapping

Every spec step is named `Scenario N Step M: ...` after the step below it performs, so the
two files read side by side. Two prose steps have no step of their own, and both are
covered — this is written down so a reader counting numbers does not take the gaps for
missing coverage:

- **Scenario 1 Step 1** (confirm no lanes column is set) is the entry state, asserted in
  the spec's `Setup` step and again by Scenario 1 Step 2's single-lane baseline.
- **Scenario 4 Step 5** (delete the probe layout even on failure) is teardown, so it runs
  in the spec's `finally` rather than as a numbered step — which is the point: a named
  step would be skipped on a red run, and this must not be.

Scenario 4 has no Step 5 in the spec's numbering for that reason; the sequence there is
1, 2, 3, 4, 6.

## Automation notes

- PRODUCT FACTS ARE NOT REPEATED HERE. Selectors, slugs, caps, defaults, menu
  inventory, drag rules, formatting and settle events live in the refdoc
  .claude/skills/grok-browser/references/viewers/tile.md; this section carries only
  test decisions and cites refdoc SECTION NAMES.
- TARGET LAYER playwright: the lane DRAGS need real input — a dispatched event does not
  reach that path, and it fails the same way the allowDragBetweenLanes=false negative case
  does, so a substitution there would look like a pass (refdoc: "Drag between lanes",
  "Selection and current-row rendering"). The property-grid toggles, the View | Layout |
  Save to Gallery leaf and the ribbon Save are driven for real for a different reason: they
  are the surfaces under test, and an apitest round-trip would not exercise them at all.
- DRIVEN FOR REAL: the toggle TRANSITIONS UNDER TEST (Scenario 1 Step 6 off→on,
  Scenario 2 Step 4 off / Step 5 on) as property-grid checkboxes
  (input[name="prop-view-show-empty-lanes"] / prop-view-allow-drag-between-lanes, via
  ensurePropertyCategory + setPropertyGridCheckbox); the lane drags as page.mouse
  gestures; the layout SAVE as the View | Layout | Save to Gallery top-menu leaf
  (driveTopMenuLeaf); the project SAVE as the ribbon Save (saveProjectViaUI).
- The one-time enablement of Allow Drag Between Lanes before Scenario 2's first drag is
  a props PRECONDITION, not a gesture: the default is not assumed, and the assert is
  that the drag WRITES a cell.
- SUBSTITUTED ACTUATION — each entry names why the surface is unavailable and which
  other channel carries the assert (a spec comment is read by no gate; only this
  section counts):
  - explicit lane LIST (props.lanes): no UI surface at all (refdoc: "The explicit lane
    list (`lanes`) has no UI surface"). ASSERT = lane DOM order, so not a prop echo.
  - lanes COLUMN (props.lanesColumnName): the only surface cannot address one named
    column (refdoc: the `prop-lanes` row under "Property grid (F4 / Context Panel)",
    and "Recon limitations / still open"). ASSERT = lane DOM structure.
  - RACE categorical FILTER (applyCategoricalFilter): a REGISTERED helper is still a
    substitution — it works through fg.updateOrAdd and clicks no filter checkbox, so
    it must never be read as driving the panel. ASSERT = tiles surviving per lane,
    cross-checked against df.filter.trueCount. KNOWN GAP: the panel's own categorical
    checkbox stays unexercised in this section.
  - grid CELL EDIT (df.set on AGE): no DOM cell editor, and the write lands on the
    same propagation path as a UI edit (refdoc: "Cell edits reach the tile without
    re-entering Edit Form").
  - column RENAME (df.columns.byName().name): the header rename is a canvas gesture
    and the assignment is path-equivalent (refdoc: "Column rename follows the DOM").
  - calculated-column ADD (df.columns.addNewCalculated) and the two cap-fitting
    REMOVEs (df.columns.remove): fixture shaping only, seating the frame on the field
    cap (refdoc: the cap and the last-added rule under "Tile internals (the read
    channel)"). ASSERT = the promoted field's value vs the grid cell text.
  - layout APPLY (tv.loadLayout): no captured selector for one gallery layout, so the
    apply leg uses the JS-API path; the SAVE leg it guards is the real menu.
  - project REOPEN (grok.dapi.projects.find + open): no captured selector for a
    per-project Browse node; the JS-API open reaches the same project-open path. The
    project SAVE it depends on is real.
- FIXTURE — promoting a newly added column onto the card needs TWO freed slots; freeing
  one leaves COMPUTED_H off the card and the step vacuous (refdoc: "Tile internals
  (the read channel)" for why).
- WITNESS DESIGN — tile reads go through the tile's own inputs, and every value is
  compared with the grid's displayed cell text, never with the raw cell value
  (refdoc: "Tile internals (the read channel)").
- WITNESS DESIGN — tile counts: assert per-lane presence/absence of a specific
  row's tile, never a total, unless the row set fits the viewport and the step
  states that precondition (refdoc: "Lanes rebuild and virtualization").
- WITNESS DESIGN — settle: gate every tile read on the viewer's own render event,
  never on a fixed sleep (refdoc: "Lanes rebuild and virtualization").
- GUARD — assert non-empty on BOTH ends before comparing two read values: '' === ''
  passes vacuously.
- GUARD — a row lookup that finds nothing leaves its -1 / '' seed, and a seed satisfies
  equality and inequality alike; guard the seed out BEFORE the comparison using it.
- LANE HEADER TEXT IS ASSERTABLE IN THIS FIXTURE: the lanes column used here (RACE)
  carries no cell renderer, so the renderer caveat in refdoc "Lanes rebuild and
  virtualization" does not reach these steps.
- DEDUP: property round-trips (lanesColumnName, lanes, tilesFont, etc.) are
  covered by public/packages/ApiTests/src/ai/viewers/tile-viewer-js-api.ts and
  are NOT repeated here. Generic Row Source semantics belong to a shared viewer
  scenario.
- TEARDOWN: delete probe layouts and projects in finally, red runs included; assign the
  project NAME BEFORE the save, or the finally block has no handle to delete by.
