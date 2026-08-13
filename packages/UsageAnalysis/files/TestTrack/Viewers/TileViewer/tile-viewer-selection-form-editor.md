---
feature: tileviewer
realizes_atlas:
  - tileviewer.cp.selection-classes-and-form-editor
  - tileviewer.int.show-selected-rows-is-viewer-local
  - tileviewer.int.autogenerate-is-a-state-flag
realizes:
  - viewers.tile-viewer
priority: p1
target_layer: playwright
coverage_type: regression
produced_from: atlas-driven
related_bugs:
  - id: GROK-19950
    status: fixed
  - id: GROK-19983
    status: fixed
  - id: GROK-19016
    status: fixed
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
expected_results:
  - anchor: "Scenario 1 Step 2"
    expectation: >-
      After a plain click on tile for row 2, that tile carries the class
      d4-current; df.currentRowIdx equals 2; df.selection.trueCount is unchanged
      (zero). No other tile gains d4-current.
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      After Ctrl-click on the tile for row 4, that tile carries d4-selected;
      df.selection.get(4) is true. After Ctrl-click on the same tile again, the
      d4-selected class is removed from it and df.selection.get(4) is false. The
      toggle is confirmed in both directions.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      After Shift-click on the tile for row 6, tile for row 6 carries
      d4-selected and df.selection.get(6) is true. Tiles for rows BETWEEN the
      previously selected tile and row 6 do NOT carry d4-selected and their
      df.selection entries are false — Shift-click is additive for that single
      row, not a range.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After Ctrl+Shift-click on the tile for row 6, the d4-selected class is
      removed from that tile and df.selection.get(6) is false. The previously
      Ctrl-selected tile still carries d4-selected and df.selection agrees.
  - anchor: "Scenario 2 Step 2"
    expectation: >-
      After setting showSelectedRows=false, the lanes host
      (.d4-tile-viewer-lanes-host) gains the class d4-tile-viewer-hide-selected;
      the computed background-color of a tile carrying d4-selected equals the
      computed background-color of an unselected tile (the selected highlight is
      neutralised); df.selection.trueCount is UNCHANGED; every selected tile
      still carries the class d4-selected.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      After setting showSelectedRows=true again, the
      d4-tile-viewer-hide-selected class is removed from the lanes host and the
      computed background-color of a d4-selected tile differs from that of an
      unselected tile (the selected highlight is restored).
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      After setting rowSource=Selected, the lanes host carries
      d4-tile-viewer-hide-selected regardless of the showSelectedRows toggle;
      the tiles rendered are only those whose df.selection is true. Setting
      rowSource back to Filtered removes the class again, and df.selection.trueCount
      is unchanged across both switches — the suppression lives on the viewer, not
      on the data. The computed opacity of the Show Selected Rows property cell
      goes 1 → 0.5 → 1 across the same round trip (the dependsOn dims it) while
      the control remains operable; no disabled attribute ever appears.
  - anchor: "Scenario 3 Step 2"
    expectation: >-
      The Edit Form designer opens as a separate view; sketchState['table']
      equals the demog frame's name, confirming the designer opened on the
      correct table.
  - anchor: "Scenario 3 Step 2b"
    expectation: >-
      The designer's EDIT button opens the Select columns... dialog; the dialog
      carries a Search field, the All and None bulk toggles, and a checked
      counter, and that counter equals the number of value fields the card
      currently renders (the 10-field cap over demog's 11 columns). Cancel closes
      the dialog leaving the column set unchanged and the designer still open.
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      After removing one field and clicking CLOSE AND APPLY, every rendered tile
      has lost the removed field's value input (input[name="input-<COLUMN>"] is
      absent for that column); the remaining field value inputs are still
      present and each equals the column's display string for that row;
      autoGenerate reads false AND sketchState.formDesigned reads true — the
      viewer has transitioned to the designed state on both channels. The
      console is clean.
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      Read the designer's host names as two sets — VALUE hosts (those holding
      input[name^="input-"]) and LABEL hosts (those holding
      input.d4-sketch-column-name) — because the channels move independently and a
      single total count cannot say which one changed. After clicking a label host
      with a real Playwright click and pressing Delete, the label set has lost
      exactly that name and the value set is unchanged. After clicking RESET, BOTH
      sets equal what they were when the designer opened: the label deleted in this
      session is back, while any value removed and APPLIED by an earlier session
      stays absent — RESET reverts the unapplied edits of the current session, it
      does NOT rebuild the default form. After CLOSE AND APPLY the reverted label
      renders on every tile, that column's value input is still present and equals
      the column's display string, and the viewer is still designed on both
      channels (autoGenerate false, sketchState.formDesigned true). The console is
      clean throughout (GROK-19016 regression guard).
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      After layout save and re-apply, the tiles show the designed field set: the
      value host removed in this step — chosen from the value hosts actually
      present, never a hardcoded column, since earlier steps reshape the set — is
      still absent from every tile. The console is clean.
  - anchor: "Scenario 4 Step 2"
    expectation: >-
      In the auto-generated state (freshly added viewer), right-click a tile field's
      value input and pick the column-removing leaf from the field menu — on this
      build that leaf is "Remove" ([name="div-Remove"]); there is no "Delete" leaf
      and no confirmation dialog. The deleted column is absent from the dataframe
      column list. Every tile has lost that column's field (its value input and
      label input are absent). The remaining tile fields are still present and
      equal their columns' display strings. autoGenerate reads true (the form
      regenerated — no user work was lost).
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      In the designed state (after CLOSE AND APPLY in Scenario 3 or a fresh Edit
      Form apply), right-click the value input of a tile field that is still visible
      and pick the same column-removing leaf, "Remove" ([name="div-Remove"]).
      The column deleted must be one the tile actually renders at that moment (read
      it off the current field set), the field set must be non-empty, and the frame
      must hold columns that have NO field — otherwise there is no freed slot and
      nothing to contrast. That column is again removed from the dataframe, and its
      field disappears from every tile: that happens in BOTH states, so it is not
      the discriminator. The mark of the designed state is the ABSENCE OF REFILL of
      the freed slot — not the survival of the deleted column's field. Concretely:
      no column that had no field before the delete gains one, and the surviving
      fields keep both their composition and their DOM order. autoGenerate reads
      false and sketchState.formDesigned reads true.
realized_as:
  - tile-viewer-selection-form-editor-spec.ts
---

# Tile Viewer — Row-state rendering and Edit Form designer

## Setup

1. Close all open views and tables.
2. Open System:DemoFiles/demog.csv (11 columns: USUBJID, AGE, SEX, RACE, DIS_POP,
   HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY; RACE has 4 categories; SEX has
   2; CONTROL is boolean). The auto-generated card is capped at 10 fields and fits 10
   of the 11, leaving SEVERITY without a field. Scenario 4 does NOT rely on that: it runs
   on its own fixture frame and grows it past the cap to manufacture the uncovered set.
3. Add a Tile Viewer to the current table view.
4. Confirm the viewer is in the auto-generated state: tiles are rendered for every
   row using the default field layout and no custom form design has been applied.
   This is the entry condition for Scenario 4.

## Scenarios

### Scenario 1: Click semantics and modifier-key selection on tiles

Row state is rendered through two independent channels — the tile's visual
selection highlight and the underlying dataset selection. Every step cross-checks
both channels.

Steps:
1. Note the visible tiles in the lane and identify tiles corresponding to rows 2,
   4, and 6 (or substitute the first three visible row indices if the display
   differs). Confirm that initially no tile is highlighted as current or selected
   and the dataset has zero selected rows.
2. Click (plain left-click) on the tile bound to row 2. Verify the click result.
3. Hold Ctrl and click the tile bound to row 4 (Ctrl-click adds row 4 to the
   selection). Verify row 4 is selected. Then Ctrl-click the same tile again.
   Verify the toggle removes the selection. (Assert both directions.)
4. Ctrl-click row 4 to re-select it. Then Shift-click the tile bound to row 6.
   Verify row 6's selection and that no intervening rows are selected.
5. Ctrl+Shift-click the tile bound to row 6. Verify row 6 is cleared from the
   selection.

### Scenario 2: showSelectedRows toggle and rowSource forced suppression

Selection presentation is controlled by the **Show Selected Rows** property and,
independently, by the **Row Source** property. The dataset selection is never
changed by either. Two independent triggers are exercised.

Steps:
1. Select two tiles using Ctrl-click (rows 4 and 6 from Scenario 1, or any two
   visible rows). Confirm both tiles are visually highlighted as selected and the
   dataset reports 2 selected rows.
2. Set the Tile Viewer's **Show Selected Rows** property to off. Verify the selection
   presentation change.
3. Set **Show Selected Rows** back on. Verify the highlight is restored.
4. Set the Tile Viewer's **Row Source** property to **Selected**. The viewer now
   shows only the two selected rows as tiles.
5. Inspect the Show Selected Rows control in the property panel and verify the
   forced suppression and the property cell opacity.

Teardown: set Row Source back to All (or Filtered) after this scenario.

### Scenario 3: Edit Form designer — field removal, label deletion, RESET, and layout round-trip

The Edit Form designer is the only surface that configures a tile's fields. It is
opened from the viewer's context menu — reached by focusing a tile and pressing the
context-menu key, because a mouse right-click lands on the field menu, which has no
**Edit Form...** entry — and it opens a separate view, not a dialog. Selecting a
field inside the designer requires a real mouse click: a synthetic event does not
select it, so a subsequent Delete key press has no effect.

Steps:
1. Focus a tile and press the context-menu key, then choose **Edit Form...** from
   the viewer menu. The designer opens as a separate view.
2. Confirm the designer opened on the demog table and that the viewer is currently
   in auto-generated mode (no custom form design has been applied yet).
2b. With the designer still open, click its **EDIT** button to open the **Select
   columns...** dialog. Confirm the dialog opens with a Search field, the **All** and
   **None** bulk toggles, and a checked counter, and that the counter equals the number
   of value fields the card currently renders (the 10-field cap over demog's 11 columns).
   The per-column checkboxes are canvas and are not touched. Click **Cancel** to close the
   dialog without changing the column set, leaving the designer open for the next step.
3. Click one field in the designer to select it (the field gains a selection
   highlight) and note the column name, then press the **Delete** key to remove it
   from the form. Click **CLOSE AND APPLY** and wait for the viewer to re-render.
   Verify the field removal and the state transition.
4. Open the viewer's context menu again (focus a tile, press the context-menu key)
   and choose **Edit Form...**. Note the field set the designer opens with — it is
   the applied state, not the factory default. Click the label area of one of the
   remaining columns to select it, then press the Delete key to remove the label. Click
   **RESET** and verify the field set returns to what the designer opened with —
   the label deleted a moment ago comes back, but the field removed and applied in
   Step 3 does not. Then click **CLOSE AND APPLY** and verify the reverted state is
   what the tiles render.
5. Save the current view layout and keep a handle on it. Then close the Tile Viewer
   and add a Grid viewer, so re-applying the layout has to restore the designed Tile
   Viewer and drop the foreign Grid. Re-apply the saved layout to the view and wait
   for re-render. Verify the designed field set is preserved. Delete the probe layout
   in teardown.

### Scenario 4: GROK-19983 two-state pin — same gesture, opposite outcomes

The same right-click > Remove gesture on a tile field produces two documented
outcomes depending on which state the viewer is in. Both must be exercised in the
same run, through the real menu gesture — that gesture is what the ticket is about.

Steps:
1. Start with a freshly added Tile Viewer (or close and re-add to reset state).
   Confirm the auto-generated entry condition: tiles are rendered using the default
   field layout and no custom form design has been applied.
2. Right-click the value of a tile field that is currently rendered and choose
   **Remove** from the field menu — this deletes the underlying column from the
   dataset (the menu has no **Delete** entry on this build, and no confirmation
   dialog appears). The column leaves the dataset, its field disappears from every
   tile, and the remaining fields still show their columns' values: the form
   regenerated, so no user work was lost. (On a card that is full to its field cap
   the freed slot is also refilled from a column that had none — that half of the
   contrast is proved in tile-viewer.md, which runs on a wider frame.)
3. Now open Edit Form, remove one additional field, click **CLOSE AND APPLY** to
   enter the designed state. Confirm the viewer is now in designed mode (the custom
   form has been applied and no further auto-generation will occur).
4. Right-click a tile field that is still visible (one that was NOT removed in Step 3)
   and choose **Remove**. The column leaves the dataset and its field disappears from
   every tile here too — that happens in both states. What marks the designed state is
   that the freed slot stays empty: no column that had no field gains one, and the
   surviving fields keep their composition and their order on the card.

Teardown: restore any deleted columns by closing and reopening demog, or by
undoing through the dataframe history if available.

## Step-to-spec mapping

Every spec step is named `Scenario N Step M: ...` after the step below it performs. Four
prose steps have no step of their own because they are preparatory or observational and are
performed inside the adjacent named step; written down so a reader counting numbers does
not take the gaps for missing coverage:

- **Scenario 2 Step 4** (set Row Source to Selected) is performed inside `Scenario 2
  Step 5`, which asserts the forced suppression the setting causes.
- **Scenario 3 Step 1** (focus a tile, press the context-menu key, choose Edit Form...) is
  performed inside `Scenario 3 Step 2`, which asserts the designer opened on the right table.
- **Scenario 4 Step 1** (start from a freshly added viewer in the auto-generated state) is
  the entry state of the SCENARIO 4 FIXTURE, guarded inside `Scenario 4 Step 2` on both
  channels — not by the spec's `Setup` step, which reads the demog viewer.
- **Scenario 4 Step 3** (enter the designed state via Edit Form → CLOSE AND APPLY) is the
  transition `Scenario 4 Step 4` needs, so it happens inside that step, immediately before
  the freed-slot assert that depends on it.

One spec step has no prose step number of its own: `Scenario 3 Step 2b` was inserted while
the designer is already open, and its anchor carries the same name.

## Automation notes

- PRODUCT FACTS ARE NOT REPEATED HERE. Selectors, slugs, designer button names, menu
  inventories, state-model reads, caps and defaults live in the refdoc
  .claude/skills/grok-browser/references/viewers/tile.md; this section carries only
  test decisions and cites refdoc SECTION NAMES.
- ACTUATION CHANNEL PER STEP. This scenario is an integration-layer link, not the
  section's UI-smoke carrier (the monolith carries UI smoke), so it does not declare
  pyramid_layer: ui-smoke and its body must promise only channels it actually drives:
  - The step-5 opacity read is a DIFFERENT signal, never a substitute for pressing the
    checkbox in step 2 — the "panel already read elsewhere" argument does not transfer.
  - Scenario 2 step 2 (Show Selected Rows off) — DRIVEN FOR REAL through the property-grid
    checkbox (input[name="prop-view-show-selected-rows"], class-1 per refdoc "Property
    grid (F4 / Context Panel)"), reached with openViewerGear + ensurePropertyCategory +
    setPropertyGridCheckbox. ASSERT = the product consequence (suppression class on the
    lanes host, neutralised selected-tile background, unchanged DataFrame selection),
    never the property value.
  - Scenario 2 step 3 (Show Selected Rows back on) — the same real checkbox back on; the
    step proves the presentation change reverses.
  - Scenario 2 step 4 (Row Source = Selected) — SUBSTITUTED by a property set although a
    real surface exists (refdoc: "A choice-valued row is a label until its cell is
    clicked"): the subject is the suppression Row Source forces on Show Selected Rows,
    and step 5 already reads the panel's reaction (the dimmed cell) to that change, so
    driving the panel adds no separate signal.
  - Scenario 3 step 5 (layout round-trip) — SUBSTITUTED by the layouts API (save, corrupt
    the view, re-apply, delete): the subject is whether the FORM COMPOSITION survives the
    cycle; the Save Layout button and the Layouts gallery are generic application chrome,
    out of this viewer's scope, and would add gallery-navigation flake to an assert about
    field sets.
- TARGET LAYER playwright: the tile click semantics and the designer's field selection
  need real mouse and keyboard input (refdoc: "Clicks need REAL input" under "Selection
  and current-row rendering", and "Editing fields on the sketch canvas — class-1"); the
  interactions are multi-step UI flows, not a JS API round-trip.
- SCENARIO 1 CHANNEL: real Playwright mouse clicks with real modifier keys, never
  dispatched events.
- WITNESS DESIGN — tile reads go through the tile's own inputs (refdoc: "Tile internals
  (the read channel)").
- SPEC MUST KEEP, NOT YET DONE — a boolean field must be compared on `.checked` against the
  raw cell value, never on `.value` against the display string. Scenario 3 Step 4 currently
  does the latter on a target that resolves to the boolean column, and it passes: so either
  the bool input mirrors the display string into `.value` — undocumented, and the refdoc says
  to read `.checked` — or the target is not the bool column. Unverified either way; do not
  "fix" the compare or the refdoc on inference. Failure direction here is red, not green.
- WITNESS DESIGN — settle: poll for the condition the step expects rather than accepting the
  first render event, and never assert on a lagging read. Where a fixed wait remains, the step's
  positive conjunct still falsifies a stale read, so the failure direction is red, not green
  (refdoc: "Lanes rebuild and virtualization").
- GUARD — before comparing the rendered tile set with the selected rows under
  rowSource=Selected, assert the selection is non-empty: zero tiles versus zero selected
  rows passes on a broken viewer.
- SELECT COLUMNS DIALOG — the step opens it from the designer, asserts the checked
  counter against the rendered field count, and leaves through CANCEL: never OK and
  never a per-column toggle, because either one reshapes the composition the later
  ladder steps assert on (refdoc: "`Select columns...` dialog (the column list)").
- COLUMN-REMOVING GESTURE IS REAL: the step right-clicks the field's value input and
  clicks the removal leaf named in refdoc "Field (column) menu — what a real
  right-click opens"; a programmatic removal would reach the same listener while
  leaving the control unexercised. Never assert a leaf COUNT, only the leaf needed.
- VIEWER-MENU OPENER: the Edit Form steps focus the tile's sketch and press the
  ContextMenu key; a right-click opens the FIELD menu instead, so a wait for the
  viewer popup after a right-click never returns (refdoc: "Context Menu" and its
  right-click pitfall).
- WITNESS DESIGN — the dependsOn signal on Show Selected Rows is the property cell's
  computed opacity, never a disabled attribute and never a value flip (refdoc:
  "`dependsOn` dims, does not disable").
- PROPERTY-PANEL PROCEDURE: expand a collapsed category with the four-argument
  ensurePropertyCategory before touching its row, probe on the category header row,
  and read a row's style with getComputedStyle inside evaluate rather than behind a
  visibility gate (reasons: refdoc "Property grid (F4 / Context Panel)").
- STATE TRANSITIONS BY STEP: Scenario 3 step 3 and Scenario 4 step 3 each move the
  viewer from auto-generated to designed, and both channels are asserted at those
  points (refdoc: "Auto-generated ↔ designed state model").
- SCENARIO 3 FIELD SELECTION: select a designer host with a PLAIN real Playwright
  click — no modifier — before pressing Delete; why a synthetic event is not enough is
  in refdoc "Editing fields on the sketch canvas — class-1". Step 3 sequence: click
  the VALUE host, press Delete, click CLOSE AND APPLY, then gate on the removed field
  leaving the tiles. Step 4 acts on a LABEL host instead, so the two steps exercise
  the two independent channels.
- WITNESS DESIGN — the clean-console guards (GROK-18230, GROK-19016) fail on ANY console
  error except one excluded BY NAME: Chrome's "Permissions policy violation:
  compute-pressure" warning, browser policy noise that lands in an arbitrary step's
  window. Exclude that one string; narrowing the guards to the crash signature instead
  would hide a regression behind the same green.
- WITNESS DESIGN — RESET is asserted on the two host-name SETS split by channel (value
  hosts vs label hosts), never on one summed host count, which cannot say which
  channel came back (refdoc: "Editing fields on the sketch canvas — class-1").
- WITNESS DESIGN — the designed-state discriminator is the ABSENCE OF REFILL of the
  freed slot, never the survival of the deleted column's field, and never "the field
  set is unchanged" (refdoc: "Auto-generated ↔ designed state model"). The no-refill
  assert is vacuous unless the frame carries more columns than the field cap and the
  excluded set is proved non-empty first, so this scenario grows its fixture past the
  cap before the contrast step.
- DESIGNER TARGETS ARE DYNAMIC: pick the host to act on from the field set as it stands
  at that moment, never a hardcoded column. Scenario 3 is a ladder (Step 3's removal
  survives Step 4's RESET) and Scenario 4 Step 4 reshapes the fixture before the column
  delete, so a fixed name addresses a host that is gone and the locator hangs.
- PROP ECHO FORBIDDEN: never set a property and assert only that the same property
  reads back. Every assertion names a product-computed consequence (DOM class,
  computed style, DataFrame state, lane structure, or clean console).
- DEDUP: property round-trips (lanesColumnName, tilesFont, etc.) are covered by
  public/packages/ApiTests/src/ai/viewers/tile-viewer-js-api.ts and are NOT
  repeated here. Generic Row Source semantics beyond the showSelectedRows
  interaction belong to a shared viewer scenario.
