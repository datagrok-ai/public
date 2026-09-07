---
feature: trellisplot
realizes_atlas:
  - trellisplot.cp.click-to-filter
  - trellisplot.int.click-cell-filter-select-events
  - trellisplot.int.onclick-filter-panel-collab
realizes:
  - viewers.trellis-plot
priority: p1
target_layer: playwright
boot_lane: local
coverage_type: regression
related_bugs:
  - id: GROK-17708
    status: fixed
  - id: GROK-17711
    status: fixed
  - id: GROK-18144
    status: fixed
  - id: GROK-17714
    status: fixed
  - id: GROK-19790
    status: fixed
  - id: GROK-13205
    status: fixed
realized_as:
  - trellis-plot-click-to-filter-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      df.filter.trueCount drops to exactly the row count matching the clicked
      cell's SEX and RACE category combination (strictly below the full demog
      row count); dataFrame.rows.filters contains one 'SEX: <value>' and one
      'RACE: <value>' entry.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      d4-trellis-plot-inner-viewer-clicked fires on mouse-down with a catIndexes
      payload; d4-trellis-plot-current-cell-changed fires with
      args.matchCondition identifying the clicked cell's category combination.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      Inner cell canvases are unchanged from the pre-click state (no visual
      alteration from the On Click mode switch alone — GROK-17708 guard).
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      The panel-only value is measured with NO trellis cell active — the cell
      clicked in Step 5 is dropped with Escape first, and every row is visible
      again before the panel filter is added (resetting the table's filter from
      outside the trellis does not drop the cell: the trellis re-applies it on
      the next recompute, and the panel-only reading would silently be an AND).
      Then df.filter.trueCount after the cell click equals (panel-only
      trueCount) AND (cell's row count), which is strictly below both the
      panel-only value and the cell-only value; the filter is the AND of both
      constraints.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      ESC on the trellis plot removes only the trellis contribution;
      df.filter.trueCount returns to exactly the panel-only value while the
      Filter Panel filter remains active and dataFrame.rows.filters no longer
      contains the trellis 'COL: category' entries. Immediately before the
      Escape the count is strictly below the panel-only value, so there is a
      trellis contribution to remove and the step is not graded on an already
      unconstrained frame.
  - anchor: "Scenario 1 Step 11"
    expectation: >-
      Changing the X split column drops the previous trellis contribution;
      df.filter.trueCount returns to the panel-only value (same as after ESC).
      Immediately before the split column changes the count is strictly below
      the panel-only value, so a contribution was demonstrably there to drop.
  - anchor: "Scenario 1 Step 12"
    expectation: >-
      Removing the Filter Panel filter through the panel itself restores
      df.filter.trueCount to the full demog row count (it was strictly below it
      immediately before), and no 'DIS_POP: <value>' entry survives in
      dataFrame.rows.filters.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Setting On Click = Filter by driving the Context Panel's 'On Click'
      property row — a different actuation path from the property assignment
      Scenario 1 uses — works identically. Row Source is parked on 'Filtered'
      through the same panel first and reads back 'Filtered' before the On Click
      commit, so the auto-correct has something to correct: Scenario 1 already
      moved it off 'Filtered', and an assert that it is "not Filtered" against a
      frame where it already is not could not fail. GROK-17711 guard: the panel
      displays 'Filter', rowSource moves from 'Filtered' to 'All', and a
      subsequent corner-click on a cell drops df.filter.trueCount to the cell's
      row count and registers 'COL: category' entries in dataFrame.rows.filters.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      Entry state of the ladder, established rather than assumed: On Click is
      Filter, no rows are filtered out and no 'SEX:'/'RACE:' entry is in
      dataFrame.rows.filters. The cell clicked in Step 4 must be cleared through
      a split-column change — neither Escape nor On Click = None clears it, and
      returning to Filter re-applies its constraint with no click at all. Then,
      walking Row Source through every value the Context Panel offers: each
      value commits and reads back, the grid keeps all 8 cells (Pack Categories
      off), no value leaves a trellis filter contribution behind
      (df.filter.trueCount stays at the full row count and rows.filters carries
      no 'SEX:'/'RACE:' entry), the Filter + Filtered combination is never left
      standing (choosing Row Source = Filtered flips On Click to None), and no
      page or console error is raised across the ladder. That flip is pinned on
      the rung IMMEDIATELY BEFORE the Filtered one, whose On Click must still
      read Filter: the ladder's entry state is several rungs away and nothing in
      between constrains On Click, so a 'None' reading on the Filtered rung
      would otherwise be as well explained by an earlier rung having quietly
      dropped Filter as by a live correction. On top of that floor, the
      GROK-13205 symptom ITSELF — cells showing IDENTICAL or EMPTY data after a
      row-source change: for every row source whose rows span more than one
      category combination (the span is computed from the live frame: df.filter,
      df.selection, df.currentRowIdx, df.mouseOverRowIdx), the per-cell canvas
      hashes are not all equal AND no cell that still receives rows is blank (no
      canvas child, or a canvas of one flat colour). Row sources whose rows span
      one combination or none (CurrentRow, MouseOver*, and any source the
      current frame state leaves narrow) are exempt — degeneracy is the correct
      rendering there. Driven-guard: the per-cell hash signature differs between
      at least two rungs, so the assert cannot pass on a grid the row source
      never reached.
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      df.selection.trueCount equals the clicked cell's structural row count; the
      +/- viewport configuration (category counts) is unchanged — GROK-17714
      guard.
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      Ctrl+click on a filter-emptied cell (structural rows exist, all hidden by
      the Filter Panel filter on SEX): df.selection.trueCount grows by exactly
      the cell's structural row count — the selection is dataset-level and
      ignores the filter. The cell must still be in the grid when it is clicked
      (all 8 cells present). The SEX filter is ADDED to the panel first and the
      panel is opened afterwards: Scenario 1 Step 12 emptied the filters group,
      and an emptied group is not repopulated on its own.
  - anchor: "Scenario 3 Step 7"
    expectation: >-
      Ctrl+click on a structurally empty cell (Pack Categories off;
      Asian|Critical combination has no dataset rows): df.selection remains
      unchanged — Ctrl is additive and an empty cell contributes no rows
      (GROK-19790 regression guard). Graded together with a driven-guard that
      the click ARRIVED: afterwards exactly ONE cell carries the current-cell
      highlight and it is the empty cell, which is a different cell from the one
      that seeded the selection. Without it an unchanged selection is also what
      a click that never reached the trellis produces. The highlight is the
      witness rather than the current-cell-changed payload: the payload is built
      from the cell's category indexes — the same quantity the regression path
      degenerates — while the highlight moves before the Ctrl modifier is
      consulted at all.
  - anchor: "Scenario 3 Step 8"
    expectation: >-
      A plain click (no Ctrl) on a different populated cell replaces the
      selection; df.selection.trueCount equals that cell's structural row count
      only. The clicked cell is a genuinely different one: its grid position
      differs from the cell that seeded the standing selection AND its row count
      differs from the standing count, so the result is distinguishable from an
      additive click (which would land on the sum of the two disjoint cells) and
      from a click that did nothing (which would leave the standing count).
  - anchor: "Scenario 4 Step 3"
    expectation: >-
      Each arrow key press fires d4-trellis-plot-current-cell-changed with a new
      args.matchCondition identifying the newly focused cell's category
      combination. Novelty is graded across ALL THREE presses, not on the first
      pair: Right → Down → Left from F|Asian visits three pairwise-different
      combinations, so the payloads must form a set of three. A comparison of
      the first two payloads alone would pass a build in which ONE arrow stopped
      moving the current cell while still announcing itself. The previously
      focused cell loses the current-cell highlight. That loss is graded as a
      MOVE of a single marker, never as its presence: the highlighted cell's
      grid position is recorded before the presses, and afterwards exactly ONE
      cell carries the current-cell highlight and it sits at a different
      position. A presence read would stay true both when the highlight never
      moved and when the previous cell kept it as well.
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      With On Click = Filter active, an arrow key moves the current cell to a
      DIFFERENT category combination and re-applies the trellis filter to it:
      df.filter.trueCount equals exactly the newly focused cell's row count (not
      merely "changed and non-zero" — a filter that was reset instead of carried
      would satisfy that).
---

# Trellis Plot — Click-to-Filter and Click-to-Select

## Setup

1. Close all open views.
2. Open **System:DemoFiles/demog.csv**.
3. Add a **Trellis plot** viewer (use the toolbox or the top menu Add viewer > Trellis plot).
4. Set the X column to **SEX** (2 categories: F, M) and the Y column to **RACE**
   (4 categories: Asian, Black, Caucasian, Other) — cell count: 2 × 4 = 8.
5. Switch the inner viewer to **Scatter plot** (the type selector in the viewer's control
   panel) so that every cell has a visible inner canvas. This is setup, not an action
   under test — see Automation notes for the channel automation uses.

## Scenarios

### Scenario 1: On Click = Filter with Filter Panel collaboration and ESC reset

**Corner-click discipline:** every actuation that uses the On Click action must click
strictly at the **corner** of a trellis cell — the center is intercepted by the inner
viewer's canvas (a scatter mark eats the click); the corner always lands on the cell
border, which is always owned by the trellis.

Steps:
1. Open the **Context Panel** for the Trellis plot and locate the **On Click** property
   under the viewer's settings.
2. Set **On Click** to **None** first, then take a settle-gated per-cell canvas
   snapshot of two different cells (this is the GROK-17708 reference state — no filter
   mode active yet).
3. Set **On Click** to **Filter**. This scenario uses the settings channel; Scenario 2
   repeats the same actuation through the Context Panel's **On Click** property row,
   and the two must stay different paths (GROK-17711).
4. Immediately take another settle-gated per-cell canvas snapshot of the same two cells.
   Verify the canvases are unchanged from Step 2 — GROK-17708 guard: the mode switch
   alone must not alter rendered cell content.
5. Corner-click the cell at the **F | Caucasian** intersection (upper-left quadrant of
   that cell's border area). Observe the number of rows passing the current filter and
   check which filters are now active.
6. Verify the two trellis events fired: confirm that a click event was emitted on
   mouse-down (with the clicked cell's category indexes) and that a current-cell-changed
   event was emitted identifying the F × Caucasian combination.
7. Confirm the canvas comparison of Step 4 once more: the two snapshots — the reference
   taken with On Click = None and the one taken right after the switch to Filter, both
   before any cell is clicked — are identical. The guard is about the mode switch alone,
   so it is deliberately read before the click, not after it: a click legitimately changes
   what the cells draw, and comparing across it would grade the wrong thing.
8. The cell clicked in Step 5 is still constraining the table, and it has to go before
   the panel is measured: click a cell and press **Escape** with the trellis focused,
   then confirm every row is visible again and no SEX/RACE filter entry is left. Only
   then add a **Filter Panel** to the view. In the **DIS_POP** filter, uncheck one
   category so that a non-trivial share of rows is filtered out. DIS_POP is deliberately
   not one of the split columns: a panel filter on SEX or RACE would be swallowed by the
   cell constraint and the AND below could hold with one contributor missing. Record the
   number of rows passing the filter after the panel filter settles (call this the
   panel-only value, strictly below the full row count).
9. With the Filter Panel filter still active, corner-click the **M | Asian** cell
   (choose a cell different from Step 5 so the filter AND is non-trivial).
   Observe the new row count passing the current filter.
10. Press **Escape** while the trellis plot viewer has keyboard focus.
    Observe the row count passing the current filter immediately after.
11. Without resetting the Filter Panel filter, change the X split column to **CONTROL**
    (a boolean categorical column, 2 values). Observe the row count passing the filter
    after the grid rebuilds.
12. Confirm the row count is still below the full row count, then remove the DIS_POP
    filter in the Filter Panel (delete it from the panel, do not merely re-check its
    categories). Verify every row is visible again (no rows filtered out).

### Scenario 2: On Click = Filter driven from the Context Panel property row (GROK-17711 actuation-path parity)

Steps:
1. Return the viewer to the default state: set On Click back to **None** and reset all
   filters so every row is visible.
2. Open the viewer's **Context Panel** through the settings gear on the viewer's title bar
   and select the **Trellis** tab (the panel also carries a tab named after the inner
   viewer type; the trellis properties are only on the Trellis tab).
3. In that same property list, first open the **Row Source** row's editor and pick
   **Filtered**, and confirm the panel reads back Filtered. This is the state the
   auto-correct below acts on: Scenario 1 already moved Row Source off Filtered, so
   without parking it back there "Row Source is no longer Filtered" would be true before
   anything was done and could never report a failure. Then open the **On Click** row's
   editor and pick **Filter**. This is the second actuation path: Scenario 1 committed the
   same setting through the settings channel, and a run that drives both through one
   channel does not test GROK-17711 at all.
4. Corner-click the **F | Black** cell. Observe the number of rows passing the filter and
   check which filters are active.
5. Walk the **Row Source** property in the same property list through every value it
   offers, one at a time. Note that the property list spells these values as the raw enum
   names (`CurrentRow`, `SelectedOrCurrent`, ...), while the viewer's context menu spells
   the same values as human captions ("Current Row", "Selected or Current") — go by what
   the property list shows. Before starting: set On Click back to **None** and reset the
   filters, uncheck **Pack Categories** (so the grid keeps a cell for every category
   combination whatever the row source is), select a broad set of rows (so the
   selection-based row sources are not empty by construction), and set On Click to
   **Filter** again so the auto-correct below has something to correct. Then clear the
   cell clicked in Step 4 by changing the X split column to another column and back —
   that is the only way it goes: Escape and On Click = None withdraw the cell's effect
   but keep the cell itself, and setting On Click back to **Filter** re-applies its
   constraint without any click. Confirm the entry state before the first rung: every row
   visible, no SEX or RACE filter entry, On Click still **Filter**.
   After each value, check that the grid still shows all 8 cells, that no rows have been
   filtered out, that the active filters list carries no SEX or RACE entry, and that no
   error is reported. Note the platform rule: picking **Filtered** while On Click is
   **Filter** switches On Click to **None** — that pair is never left standing. Read On
   Click on the value walked immediately BEFORE Filtered as well: it must still say
   **Filter** there. Only then is the switch something seen happening; read against the
   state at the start of the walk, several values back, a **None** on the Filtered value
   is just as well explained by an earlier value having quietly dropped Filter.
   Then look at what the cells actually DRAW, which is the reported defect this ladder
   guards: when the chosen row source feeds rows to more than one category combination
   (All, Filtered, Selected and the like), the cells must not all show the same picture,
   and no cell that still receives rows may be blank. A quiet ladder that raises no error
   while every cell repeats one picture is exactly the defect passing. When the chosen row
   source by its meaning feeds a single row or none (`CurrentRow`, the `MouseOver*` ones), a
   repeated or blank picture is legitimate and is NOT checked — grading it would be false
   red. Across the whole walk the picture must change at least once: if it never does, the
   row source is not reaching the drawing at all and "the cells differ" would be judging a
   grid nothing has touched.
6. Restore Row Source to **All**, Pack Categories back on, On Click to **None**, and clear
   the selection and all filters.

### Scenario 3: On Click = Select — additive Ctrl+click and the GROK-19790 empty-cell guard

Steps:
1. Set **On Click** to **None**, reset all filters (every row visible).
   Record the initial +/- viewport cell-count state (number of visible category rows
   and columns under the current SEX × RACE split).
2. Set **On Click** to **Select**.
3. Corner-click the **F | Caucasian** cell. Observe the number of selected rows shown
   in the status bar.
4. Verify the selection count and the viewport configuration (GROK-17714 guard: the
   +/- category count must be unchanged after the On Click = Select switch and the
   first click).
5. The Filter Panel is empty at this point — Scenario 1 Step 12 removed its last filter
   and the panel does not bring its column filters back on its own. Add a **SEX** filter
   to it and keep only **M**, so that the F | Caucasian cell's rows are excluded by the
   filter (but the structural combination still exists). Check
   that the grid still shows all 8 cells — the F column has to survive for the next step
   to be about a filter-emptied cell rather than about a different cell.
6. Ctrl+corner-click the **F | Caucasian** cell (which now appears filter-emptied —
   its rows are hidden by the panel filter, but the cell is still visible and clickable).
   Observe the number of selected rows shown in the status bar. Afterwards remove the SEX
   filter from the panel and confirm every row is visible again.
7. With **Pack Categories** turned **off** (set via Context Panel > Pack Categories
   unchecked so that structurally empty cells appear in the grid), set X to **RACE**
   and Y to **SEVERITY** so that the Asian | Critical combination has no dataset rows
   at all. Build up a selection on any populated cell first, then Ctrl+corner-click the
   **Asian | Critical** empty cell. Observe the selected row count. Note also which cell
   now carries the current-cell border highlight: it must be the empty cell just clicked,
   and no other cell may carry it. That is the proof the click actually reached the
   trellis — "the selected row count did not change" is equally what a click that never
   arrived would show, and the empty cell is deliberately a different cell from the one
   the selection was built on, so the highlight has somewhere to move from.
8. Corner-click (no Ctrl) another populated cell to perform a plain select — and pick it
   deliberately: it must be a cell other than the one whose rows make up the current
   selection, and it must hold a different number of rows than that selection currently
   shows. Re-clicking the same cell leaves the count exactly where it was, and then a
   click that replaced, a click that added and a click that did nothing all look alike.
   Observe the selected row count.
9. Reset the selection (click an empty area or use the toolbar) and restore the
   original split (X = SEX, Y = RACE) and Pack Categories = on.

### Scenario 4: Arrow-key navigation — current-cell-changed events and the filter following the current cell

Steps:
1. Reset all filters and selection; set On Click to **None**; X = SEX, Y = RACE.
2. Corner-click a cell to give the trellis keyboard focus. Confirm the current-cell border
   highlight is on that very cell and on no other one, and note which cell it is — the
   arrow presses below are graded by where the highlight ends up relative to here.
3. Press the **Right** arrow key. Press **Down** arrow key. Press **Left** arrow key.
   Verify that a current-cell-changed event is emitted identifying the newly focused
   cell's category combination, and that the three presses report three DIFFERENT
   combinations between them — Right, Down and Left from F | Asian pass through three
   different cells, so three different combinations is what a healthy run reports.
   Checking only that the first two reports differ lets one dead arrow through: the two
   remaining presses still move the highlight and the first pair still differs.
   Check as well that afterwards the current-cell highlight is on
   exactly one cell and that cell is not the one noted in Step 2 — the highlight has to
   have MOVED. A highlight that stayed put, or one that stayed on the old cell while a
   second cell also lit up, is a failure; simply seeing "some cell is highlighted" would
   accept both.
4. Set On Click to **Filter**, corner-click a cell so a trellis filter is active (confirm
   the row count is below the full row count), then press the **Down** arrow key once.
   Note which category combination the current cell moved to, and the number of rows
   passing the filter after the key press.
5. Reset On Click to None and clear the trellis filter so every row is visible again.

## Automation notes

- CHANNEL — Setup step 5 (inner type) and the Scenario 1 On Click switch are prop writes
  (`props.viewerType`, `props.onClick`): setup and a deliberate channel split, not the action
  under test. The control-panel type selector is the manual path.
- CHANNEL — Scenario 2 Step 3 commits both `rowSource` and `on-click` through the Context Panel
  property editor (`setTrellisChoiceViaPanel`), then asserts `rowSource === 'All'`, the exact
  landing the refdoc Interaction table pins. That is the GROK-17711 split against Scenario 1.
- CHANNEL — On Click actuations click a corner offset (e.g. `{x: 4, y: 4}`) of
  `.d4-trellis-plot-cell`, never the centroid (refdoc: pitfall 11).
- CHANNEL — the ONE exception to that rule is Scenario 1 Step 6, which needs BOTH events off one
  cell: `d4-trellis-plot-inner-viewer-clicked` rises only from a mousedown at the inner canvas
  CENTRE and `d4-trellis-plot-current-cell-changed` only from the corner click, and no single
  gesture raises both. So the step keeps TWO gestures on that cell — centre mousedown, then corner
  click; a single corner click leaves the inner payload null. The corner rule stands everywhere
  else and is not weakened by this.
- CHANNEL — Filter Panel steps use `openFilterPanel` + `applyCategoricalFilter` on `DIS_POP`
  (Scenario 1) / `SEX` (Scenario 3); the helper drives a categorical filter. Call order differs by
  situation: Scenario 1 Step 9 opens the panel first, Scenario 3 Step 6 adds the filter first
  (refdoc: Filter Panel — an emptied filters group is never repopulated), `openFilterPanel` in
  front of an add sitting out its 15 s `.d4-filter` wait. Selectors: refdoc: Filter Panel.
- CHANNEL — Row Source is driven and read on the property-grid `<select>`, the surface this spec
  reads (refdoc: Row Source: two surfaces); the menu captions are a different spelling.
- CHANNEL — the owning property-grid category is read off the live grid rather than hardcoded
  (refdoc: Inner Viewer Property Grid), a collapsed row being zero-size and unclickable.
  `propCategoryOf` matches the `property-grid-category` class and scopes the walk to the owning
  grid (refdoc: A `prop-category-` NAME is not a category header).
- CHANNEL — Pack Categories is driven as a property (`props.packCategories`), not through the
  Context Panel checkbox: set false for the Scenario 2 Step 5 ladder and again before Scenario 3
  Step 7, restored to true after each.
- WITNESS — `d4-trellis-plot-inner-viewer-clicked` and `d4-trellis-plot-current-cell-changed`
  listeners are attached via `page.evaluate` BEFORE each actuation; payloads are asserted after
  the event.
- WITNESS — page-context reads: `viewer.dataFrame.filter.trueCount`, `.selection.trueCount`, and
  `.rows.filters` (asserted for the presence/absence of `SEX: <value>` / `RACE: <value>`).
- WITNESS — settle-gated canvas diff (GROK-17708): a settle wait after the flip, then the same two
  cells re-read and compared. Anti-vacuity is ASSERTED, not implied by the wording "two distinct
  cells": all four readings must be non-null AND the two pre-flip readings must DIFFER FROM EACH
  OTHER. Two blank canvases read equal, and "unchanged" then holds whatever the flip did.
- WITNESS — the current-cell highlight is read as a POSITION: the index of
  `.d4-trellis-cell-current` among the viewer's cells, plus the count carrying the class.
  Scenario 4 Step 3 requires `length === 1` and a changed index; Scenario 3 Step 7 uses the same
  pair as its Ctrl+click delivery witness (`count === 1`, `index === emptyIdx`,
  `emptyIdx !== popIdx` asserted first). The class is the witness rather than the event payload
  (refdoc: The trellis filter contribution and the current cell are SEPARATE state): on the
  GROK-19790 path the payload's source field is the one that empties.
- WITNESS — Scenario 4 Step 3 grades payload novelty as a SET: each `matchCondition` serialized as
  `[SEX, RACE]` so key order cannot split one combination in two, with at least three distinct
  values — the number of cells Right → Down → Left visits.
- WITNESS — Scenario 1 Step 9 clears the trellis contribution by corner-click plus Escape and
  asserts `filter.trueCount === fullRowCount` with no SEX/RACE entry before touching the panel.
- WITNESS — the X-column round trip before the Scenario 2 Step 5 ladder is what clears the standing
  current cell: Escape and On Click = None withdraw the contribution but leave the cell (refdoc:
  The trellis filter contribution and the current cell are SEPARATE state).
- WITNESS — Row Source ladder (Scenario 2 Step 5): the walked values are read off the live choice
  editor and graded as an exact SET, both sides sorted since option order is not a guarantee. Pack
  Categories must be OFF (refdoc: What packing actually collapses), or a narrow source
  legitimately collapses the grid and the cell-count check means nothing.
- WITNESS — Row Source auto-correct: the ladder records `props.onClick` per rung and pins BOTH the
  `Filtered` rung at `None` and the rung before it at `Filter` (refdoc: The `onClick` ⇄ `rowSource`
  auto-correct runs in BOTH directions); if `Filtered` is offered first, the entry read stands in
  as the preceding witness.
- WITNESS — Row Source render signal (GROK-13205): the per-cell hash is Scenario 1's
  `getImageData` probe; a cell is empty when it has no `canvas` child or a canvas of one flat
  colour. Which rows a source feeds is derived live from `df.filter`, `df.selection`,
  `df.currentRowIdx`, `df.mouseOverRowIdx` — the mouse-over group is empty, the ladder never
  hovers — the direct read being unreachable from JS (refdoc: pitfall 26). Wide/narrow is not
  declared per value: a rung is graded only when its rows span at least two category combinations,
  so a source narrow through frame state is exempted for the right reason, and an offered value
  the spec cannot express fails loudly instead of dropping out silently.
- HONEST FLOOR — "blank" means one flat colour over the whole canvas, so a cell drawing chrome but
  no data still counts as painted; only the fully blank cell is caught. The identical-cells half is
  unaffected: chrome is the same in every cell.
- WITNESS — Scenario 3 Step 8 picks its target from the live frame, not by "first cell with a
  canvas": per-cell structural and selected counts come from `df.selection.getSelectedIndexes()`
  and the split columns' categories, the seeded cell is the one with a non-zero selected count,
  and the target is a populated cell with a different DOM index AND a different structural count.
  With no such cell the step fails rather than falling back on the seeded one.

### Spec must keep

- **Scenario 2 drives a DIFFERENT actuation path from Scenario 1** — the real Context Panel
  property row (`setTrellisChoiceViaPanel`), never `props.onClick`. Substituted, Scenario 2 is a
  duplicate of Scenario 1 on another cell and GROK-17711 loses its only coverage.
- Scenario 2 Step 3 parks Row Source on `Filtered` before committing On Click = Filter, through
  the panel editor — else the auto-correct assert has nothing to correct and cannot fail.
- Its post-commit assert is the exact landing `All`, not merely "not Filtered".
- The auto-correct is pinned on the rung IMMEDIATELY BEFORE `Filtered`, not on the entry state:
  keep the companion assert that the preceding rung still reads `Filter` (refdoc: The `onClick` ⇄
  `rowSource` auto-correct runs in BOTH directions), with the entry read as fallback when
  `Filtered` comes first. Alone, `filteredRung.onClick === 'None'` is as well explained by an
  earlier rung dropping `Filter`.
- Filter Panel steps use the real panel (`openFilterPanel` + `applyCategoricalFilter`), never a
  hand-rolled `df.onRowsFiltering` subscription — a private subscription closes
  `int.onclick-filter-panel-collab` against a simulation.
- A Filter Panel use following a full filter removal ADDS the filter first and opens the panel
  second (Scenario 3 Step 6; refdoc: Filter Panel) — the reverse order times out on `.d4-filter`.
  The panel stays the channel under test: `applyCategoricalFilter` goes through `fg.updateOrAdd`,
  and `openFilterPanel`'s own wait plus hover proves the panel carries it.
- Step 12 resets through the Filter Panel, never `df.filter.setAll(true)` (nor `v.resetFilters`,
  which ends in `setAll(true)`) — writing the bitset all-true and asserting it is all-true is a
  tautology no defect can fail.
- Step 12 keeps its pre-reset witness that the count was strictly below full.
- The owning property-grid category is found by the `property-grid-category` CLASS, never a
  `prop-category-` name prefix, and the walk stays inside the `.property-grid` that owns the row
  (refdoc: A `prop-category-` NAME is not a category header).
- The trellis contribution is removed the PRODUCT way before anything is measured against it —
  click + Escape in Scenario 1 Step 9, with the full-row-count and no-SEX/RACE witness BEFORE the
  panel filter is added. `setAll(true)` + `requestFilter()` does not remove it (refdoc: The trellis
  filter contribution and the current cell are SEPARATE state), so the "panel-only" baseline would
  silently be an AND.
- Scenario 2 Step 5 establishes its entry state, it does not assume it: the Step 4 cell survives
  Escape and On Click = None (refdoc: The trellis filter contribution and the current cell are
  SEPARATE state), so the ladder clears it with a split-column round trip
  (X → CONTROL → SEX) plus a `requestFilter`.
- It keeps both entry asserts (`trueCount === fullRowCount`, no SEX/RACE entry) — they are the
  reads every rung is graded on, and a surviving cell fails every rung with an unreadable cause.
- It also seeds a BROAD selection before the first rung (`df.selection.init(i => i % 2 === 0)` on an
  unfiltered frame) — else `Selected`, `FilteredSelected` and `SelectedOrCurrent` feed one
  combination or none, are classified narrow and drop out of the GROK-13205 render check silently,
  three rungs of eight losing coverage without one red.
- Scenario 3 Step 8 must click a cell the seeded selection does NOT come from — different grid
  position AND different structural row count — else the equality holds for a replace, an add and
  a no-op alike.
- It keeps all three companions: the entry asserts (candidate exists, index differs, count differs)
  and the two post-click negatives (not the previous count, not the sum of the two disjoint cells).
- **The current-cell highlight is graded by MOVEMENT of a single marker, never by presence.**
  `!!document.querySelector('.d4-trellis-cell-current')` stays true when the highlight never moved
  and when a second cell lights up beside the first. Scenario 4 Step 3 records the index before the
  arrows and afterwards requires exactly one element with the class plus a changed index — the
  render-signal index's form for `trellis_cell_click` and `trellis_keyboard`.
- **Scenario 4 Step 3 grades the payloads of ALL THREE arrows, never the first pair** — at least
  three distinct `matchCondition` values, serialized `[SEX, RACE]`. `cc[0] !== cc[1]` admits an
  arrow that stopped moving the current cell while still announcing itself.
- Scenario 3 Step 7 keeps its delivery witness beside the unchanged-selection assert — exactly one
  cell carrying `.d4-trellis-cell-current` afterwards, and it the empty cell, with
  `emptyIdx !== popIdx` asserted first; `selAfter === selBefore` is also what a click that never
  reached the trellis produces.
- Do not swap in the event payload there: on this path the class is the robust witness and the
  payload is not (refdoc: The trellis filter contribution and the current cell are SEPARATE state).
- **Scenario 4 Step 4 asserts EQUALITY with the newly focused cell's row count**, taken from the
  event's `matchCondition` — `changed && > 0` is satisfied by a filter reset rather than carried.
- Steps 10 and 11 keep their pre-actuation witnesses (`before < panel-only`) and Scenario 3 Step 6
  keeps the 8-cell check after the SEX filter — else an ESC, a split change or a click that did
  nothing still passes.
- Scenario 3 Step 6 seeds the selection on a DIFFERENT populated cell (`M | Black`) and asserts it
  non-zero BEFORE the Ctrl+click — out of an empty selection `0 + N = N`, and additive reads exactly
  like replace.
- Scenario 1 Step 9 asserts the panel-and-cell INTERSECTION non-empty before the click
  (`mAsianPanel > 0`) — zero satisfies both "strictly below" inequalities, and the AND would then be
  graded on a fully filtered frame.
- The error collector stays (`pageerror` + console errors, benign filter, window-gated per step);
  the Row Source ladder rests on it for the crash class. A global `toEqual([])` over the run would
  read ambient noise as red.
- That benign filter stays NARROW — resource / 404 / favicon / cloned-iframe ambient noise ONLY.
  `NullError: method not found` and Dart stack traces (`ProjectMeta.publish`, `project_meta.dart`)
  must COUNT as errors here whatever a sibling spec in this folder filters: they are the crash floor
  the GROK-13205 ladder rests on, and importing a "shared" filter retires it silently.
- **The Row Source ladder keeps its per-cell RENDER assert — the no-error floor does not hold
  GROK-13205 alone.** The ticket is cells showing IDENTICAL or EMPTY data, which raises no error,
  keeps all 8 cells and leaves `df.filter` untouched. Keep all five parts: (a) the fed cells do not
  all share one canvas hash; (b) no cell still receiving rows is blank; (c) both graded only where
  the rows span more than one category combination — grading `CurrentRow` / `MouseOver*` would be
  false red, their degeneracy being legitimate; (d) the ladder-level driven-guard that the hash
  signature differs between at least two rungs, without which a row source that never reaches the
  render leaves the previous frame standing and (a) passes on an untouched grid; (e) the
  classifier's own self-check — `All` must come out WIDE, since with the frame unfiltered and half
  the rows selected it cannot be anything else, and a classifier that called everything narrow
  otherwise hands (a) and (b) an empty set to grade.
- **The expected Row Source options are the PROPERTY-GRID spelling and the check stays an exact
  SET**, both sides sorted since order is not a guarantee — the two surfaces disagree on most of
  the values (refdoc: Row Source: two surfaces).
- Never weaken that to an "at least N options" floor: a build that dropped or renamed an option
  walks a shorter ladder and still passes.
- `props.rowSource` may be compared with the `<select>` (same spelling per the refdoc), but say so
  where it is done.
- **Row Source must end the ladder on `All`.** Scenario 3 filters SEX = M and still needs the F
  column present, which holds only while the trellis reads all rows; the default `Filtered` repacks
  the grid and the F|Caucasian index then addresses a different cell.
- **No literal row count in the prose** — it is read live at setup and pinned once in the spec.
