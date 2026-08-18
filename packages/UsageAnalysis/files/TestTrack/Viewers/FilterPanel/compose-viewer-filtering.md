---
feature: filters
realizes_atlas:
  - filters.cp.compose-with-viewer-filtering
  - filters.int.and-combination
realizes:
  - viewers.filters
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - compose-viewer-filtering-spec.ts
related_bugs:
  - id: github-2642
    status: fixed
  - id: GROK-18281
    status: fixed
  - id: GROK-16713
    status: fixed
scope_reductions:
  - id: SR-01
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      Scenario 3 Steps 2 and 3 build the SELECTION over the JS API:
      `df.selection.copyFrom(df.filter)` for "select the visible rows" and
      `df.selection.init(...)` for "add the second category to the selection".
      The two commands the scenario exists to test — Select > Invert and
      Select > Selection to Filter — are both driven through the real top menu
      leaf (`v.driveTopMenuLeaf`), and every claim about the panel card
      (criterion kept, canvas repainted) is read off the real UI. What is
      reduced is only how the selection those commands operate on is put in
      place.
    verdict_status: SCOPE_REDUCTION
  - id: SR-02
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      The PANEL criterion is seeded and changed over the JS API
      (`v.resetFilters` + `v.applyCategoricalFilter`) rather than by clicking the
      RACE card's category rows. This scenario's subject is the COMPOSITION of a
      panel criterion with viewer-driven filtering; the card gestures themselves
      are asserted end to end in panel-core-ladder (canvas category click) and in
      save-and-reapply-state (Step 2). The seeded criterion is verified through
      the UI at every step — the header counter reads 1 and the card reports
      exactly the seeded category.
      Scenario 1 Step 5 is the one place where the criterion CHANGE is itself the
      actuation under test rather than the scene, so both paths were measured
      there on dev 2026-08-18. With a scatter zoom already biting (354 → 222
      rows), clicking the RACE card's "Black" category-name row left the frame at
      127 filtered rows against the 157 rows that category holds in full, and the
      Scatter Plot still read `zoomAndFilter: filter by zoom` — the zoom survived
      a card-driven criterion change exactly as it survives the API-driven one,
      which is the github-2642 guard the step asserts.
    verdict_status: SCOPE_REDUCTION
  - id: SR-03
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      Viewers are added with `tv.addViewer(type)` and three of the Scenario 4
      channels are ARMED by setting a property (Histogram `filteringEnabled`,
      PC Plot `showFilters`, Bar chart `splitColumnName`) instead of through the
      property panel; the arming that carries the behaviour under test —
      On Click > Filter for the Trellis plot, the Pie chart and the Bar chart — is
      driven through the real context menu and confirmed by reading the property
      back. Every FILTERING gesture in this spec is a real mouse gesture.
    verdict_status: SCOPE_REDUCTION
  - id: SR-04
    check: E-ASSERTION-STRENGTH-01
    rationale: |
      Gesture DELIVERY has an independent probe on two of the six channels only:
      the Histogram handle is located by the canvas cursor turning `ew-resize`,
      and a Trellis cell click counts only when it leaves a
      `.d4-trellis-cell-current` cell behind. For the Pie chart click and its
      outside-click revert, the PC Plot slider double-click revert and the Trellis
      Escape revert there is no such probe: the helpers report only that the
      target was on screen and the gesture was ISSUED, and their assertion
      messages now say exactly that. Delivery for those is carried by the
      trueCount transition asserted immediately after — strictly below
      trueCount_panel for the drive, exactly trueCount_panel for the revert —
      which no undelivered gesture satisfies.
    verdict_status: SCOPE_REDUCTION
  - id: SR-05
    check: E-LAYER-COMPLIANCE-01
    rationale: |
      In Scenario 2 the layout is SAVED through the real View > Layout > Save to
      Gallery leaf, but RE-APPLIED over the JS API (`grok.dapi.layouts.find(id)` →
      `tv.loadLayout`) because the gallery's thumbnails carry nothing that
      identifies the layout this run just saved. The re-apply is awaited on
      `grok.events.onViewLayoutApplied` rather than a fixed sleep, and everything
      asserted about the result — rebuilt cards read off the DOM, restored
      criterion, restored row count, error balloons, page errors — is read from
      the real UI.
    verdict_status: SCOPE_REDUCTION
  - id: SR-06
    check: E-ASSERTION-STRENGTH-01
    rationale: |
      "The scatter-plot zoom comes back with the layout" is NOT assertable on
      this build: the layout round-trip does not restore it. Measured three times
      on dev — before the save the zoom plus the panel criterion left 49 rows;
      after Save to Gallery → loadLayout the Scatter Plot is back and its
      zoomAndFilter still reads "filter by zoom" — the spec asserts both of
      those, so a build that stops restoring the viewer, or restores it
      disarmed, fails here instead of passing on the same row count — but the
      row count settles at 72,
      the panel criterion alone, and stays there for the full 30 s the count is
      polled. The viewer's viewport is not visible in its properties either
      (xMin/xMax/yMin/yMax read null both before the save and after the re-apply),
      so there is no second channel to assert the zoom through. What the step
      asserts instead is the exact value that separates the two failure modes
      around it: trueCount EQUALS trueCount_panel after the re-apply — below it
      would mean the zoom is now restored (this reduction is then obsolete and the
      bound tightens to strictly below), above it would mean the panel criterion
      itself was lost. The pre-save count is asserted to be strictly below
      trueCount_panel so the equality is a statement about a narrowing that
      demonstrably existed. Operator ruling 2026-08-18: the scatter-plot zoom is
      not part of what a layout stores, so a round trip that comes back with the
      panel criterion alone is CORRECT and there is no defect to file. The
      equality is therefore the expected value rather than a hedge around
      unknown behaviour, and it still separates the two failure directions named
      above.
    verdict_status: SCOPE_REDUCTION
unresolved_ambiguities: []
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      Entering the block the Scatter Plot's zoomAndFilter property reads "filter
      by zoom", so the drag really routes to df.filter, and the row count on
      entry equals trueCount_panel. After the drag df.filter.trueCount drops
      below trueCount_panel and stays strictly above zero (scatter zoom
      intersects with the panel criterion via the AND-BitSet mechanism rather
      than replacing it or wiping every row); the active-filter counter in
      .d4-filter-group-header .d4-filter-indicator still reads 1 (no new panel
      card was created) and the RACE card still carries exactly the seeded
      category (the viewer did not rewrite the panel criterion).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After changing the panel criterion while zoom is still applied, trueCount
      changes to a new value that reflects both the updated panel criterion AND
      the still-active zoom (github-2642 regression guard: the new panel filter
      does not reset the viewer-driven one), stays strictly above zero, and the
      active-filter counter still reads 1. After the criterion is put back to the
      seeded category the row count is still strictly BELOW trueCount_panel — the
      zoom is untouched by either criterion change.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      Entering the block the row count equals trueCount_panel exactly (the
      Scatter Plot was closed, so the narrowing measured is the bar chart's
      alone) and the Bar Chart's onClick property reads Filter. A bar click on top of the active
      panel criterion is accepted only when it leaves EXACTLY ONE SEX category
      standing among the filtered rows, so which bar was hit is pinned rather than
      left to the run; the surviving row count then equals exactly the number of
      rows that are both the panel category AND that SEX value — the intersection
      computed from the data, not merely "some smaller number". It is below
      trueCount_panel and above zero, and the active-filter counter remains 1
      (bar-chart click does not become a panel card).
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      Immediately before the close the row count is confirmed to be strictly
      below trueCount_panel, so the zoom narrowing being released is proven to
      exist. After closing the Scatter Plot viewer, trueCount equals
      trueCount_panel exactly (the panel criterion is intact; the viewer-driven
      zoom narrowing was removed with the viewer), the active-filter counter
      reads 1 and the RACE card still carries exactly the seeded category.
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      Immediately before the close the row count is confirmed to be strictly
      below trueCount_panel, so the bar-chart narrowing being released is proven
      to exist. After closing the Bar Chart viewer, trueCount again equals
      trueCount_panel exactly (the bar-chart filtering is gone; the panel
      criterion alone remains), the active-filter counter reads 1 and the RACE
      card still carries exactly the seeded category.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      After saving the layout, closing the Filter Panel and re-applying the
      layout (with panel filter and scatter zoom active), no error balloon
      appeared at any moment inside that window and no page error was emitted
      (GROK-18281 regression guard). The save produced exactly ONE new layout
      belonging to the current user (a save that produced none or several would
      make the re-apply act on an unknown object). The closed panel is rebuilt by
      the re-apply — at least one .d4-filter card is back inside
      [name="viewer-Filters"], read off the DOM rather than through
      getFiltersGroup(), which would create a panel and manufacture its own green
      — and carries its saved RACE criterion again. trueCount after the re-apply
      is polled to EXACTLY trueCount_panel: the two bounds that matter here are
      trueCount_panel (the panel criterion alone) and the pre-save zoomed count
      (panel criterion AND scatter zoom), and it is trueCount_panel that the
      round trip lands on — a layout does not store the scatter-plot zoom, so
      restoring the panel criterion alone is the expected outcome (operator
      ruling 2026-08-18, scope_reductions SR-06). The pre-save zoomed count is
      confirmed strictly below trueCount_panel, so the equality is asserted about
      a narrowing that demonstrably existed; a count below trueCount_panel would
      mean the zoom came back through the layout, a count above it means the
      panel criterion itself was lost. The error-balloon observer is
      installed before the save and held until the restore condition is met AND
      no further balloon has appeared for a quiet interval, then disconnected —
      a balloon arriving late in the window is still seen.
  - anchor: "Scenario 4"
    expectation: >-
      Each of the four remaining viewer-driven filtering channels — the
      Histogram's own range slider, a PC Plot per-axis slider, a Trellis cell
      click and a Pie Chart segment click — behaves exactly as the Scatter and
      Bar channels do. Entering each block trueCount equals trueCount_panel (so
      the narrowing measured is that viewer's alone); after the gesture
      trueCount is strictly below trueCount_panel and strictly above zero, the
      active-filter counter still reads 1 and the RACE card still carries its
      single seeded category; after the viewer's OWN revert gesture (the
      histogram handle dragged back to the right end of its track, a PC Plot
      slider double-click, Escape, a click outside the pie) trueCount returns to
      trueCount_panel exactly with the viewer still open; and after closing the
      viewer it is trueCount_panel again with the counter at 1. Two channels prove DELIVERY
      independently: the histogram handle is located by the canvas cursor reading
      ew-resize and the block fails when no handle is found, and a trellis cell
      click counts only when it leaves a .d4-trellis-cell-current cell behind.
      For the remaining gestures (the pie click and its outside-click revert, the
      PC Plot slider double-click revert, the Trellis Escape revert) the helper
      reports only that the gesture was ISSUED onto a target that was on screen,
      and delivery is carried by the exact trueCount transition asserted right
      after it (scope_reductions SR-04). When a trellis cell click narrows nothing
      and is reverted with Escape, the revert is asserted to have put the row
      count back to the value the click started from before the next cell is
      tried, so no cell is ever measured from an unknown baseline.
  - anchor: "Scenario 3 Step 2"
    expectation: >-
      No viewer other than the Grid and the Filters panel is present, so nothing
      but the panel criterion is filtering and the row effect measured in Step 5
      is readable, and the panel-only base is re-established and confirmed equal
      to trueCount_panel on entry. Selecting the visible rows leaves
      df.selection.trueCount equal to trueCount_panel exactly.
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      The second category is real (its full-dataset count is strictly above zero)
      and disjoint from the panel category, and adding it grows the selection to
      exactly trueCount_panel + that category's full-dataset count — so neither
      the arithmetic below nor the command's effect can be satisfied by zeros.
  - anchor: "Scenario 3 Step 4"
    expectation: >-
      Driving the real Select > Invert leaf leaves df.selection.trueCount equal to
      exactly rowCount minus the pre-inversion selection size — the complement,
      not merely "a different number" — and that value is strictly above zero and
      strictly below rowCount. The filter and the selection are confirmed to
      DISAGREE at this point (df.filter.trueCount equals trueCount_panel and is
      not equal to the inverted-selection count), so the Step 5 equality cannot be
      satisfied by a command that does nothing at all.
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      After running Selection to filter, trueCount equals the inverted-selection
      row count recorded in Step 4 (the filter BitSet was replaced by a verbatim
      copy of the selection), and the panel's own category is now the one that is
      100% filtered out — the divergence between what the panel is set to and what
      the data holds, which is what makes the card check bite. The RACE card KEEPS
      its criterion — the command copies the selection onto the filter, it does not
      rewrite panel card settings — and the RACE card canvas REPAINTS: it differs
      from the frame captured immediately before the command, having been
      pixel-identical across the idle window that preceded it (GROK-16713
      regression guard).
---

# Filters — Composition of Filter Panel with Viewer-Driven Filtering

## Setup

1. Open `System:DemoFiles/demog.csv` (5,850 rows, 11 columns).
2. Open the Filter Panel by clicking the funnel icon in the ribbon toolbar.
   Wait for at least one filter card to appear in the panel before proceeding.
3. Apply a categorical filter on the RACE column: in the Filter Panel, click
   the **Asian** row in the RACE filter card. Wait for the filtered row count
   to drop below 5,850.
4. Record **trueCount_panel** — the number of rows currently passing the filter. RECORD it at run
   time from the product rather than expecting a fixed number; on demog the Asian category holds 72
   rows (the RACE categories are Caucasian 5,267, Other 354, Black 157, Asian 72). Record
   **rowCount** = 5,850.
5. Verify the active-filter counter in the Filter Panel header reads **1**.

## Scenarios

### Scenario 1: Filter Panel criterion composes with viewer-driven filtering via the AND BitSet

This scenario covers the core critical path: two filtering sources (Filter Panel + viewer zoom /
bar-chart click) must intersect (AND) without either clobbering the other.

Steps:

1. Add a **Scatter Plot** viewer to the layout.
2. In the Scatter Plot context menu (right-click the viewer or use its hamburger menu), find
   **On Zoom > Filter** and enable it so that zooming the scatter plot narrows the active filter.
3. In the Scatter Plot, drag a selection rectangle on the plot canvas to zoom into a region that
   contains fewer rows than trueCount_panel. Wait for the filtered row count to update.

4. Verify: on entry to this step, before the drag, the number of rows passing the filter equals
   trueCount_panel — the panel criterion is the only thing filtering, so whatever the drag changes
   is attributable to the zoom alone. Verify: the Scatter Plot's zoom-and-filter mode really reads
   **filter by zoom**, so the drag reaches the dataframe filter at all.

   Verify: the number of rows passing the filter is strictly less than trueCount_panel (the zoom
   narrows the already-filtered set — the two criteria intersect, they do not replace each other)
   and strictly greater than zero (the zoom narrowed the set, it did not wipe it).
   Verify: the active-filter counter in the Filter Panel header still reads **1** (the viewer-driven
   filtering did not create a new panel card and did not reset the panel criterion).
   Verify: the RACE filter card still carries exactly the seeded category (**Asian**) — the zoom did
   not rewrite the panel criterion.

5. With the scatter-plot zoom still applied, change the RACE panel criterion: in the Filter Panel,
   click the **Caucasian** row in the RACE filter card to replace Asian with Caucasian. Record
   **trueCount_after_change** — the new filtered row count.

   Verify: the filtered row count moved to a new value (trueCount_after_change differs from the
   value recorded in Step 3) AND the zoom-driven narrowing is still in effect
   (trueCount_after_change is less than the total number of Caucasian rows in the full dataset,
   meaning the zoom is still active). This is the **github-2642 regression guard**: applying a new
   panel filter must not silently reset the viewer-driven zoom.
   Verify: trueCount_after_change is strictly greater than zero, and the active-filter counter in
   the Filter Panel header still reads **1** (changing the criterion did not add a second card).

6. Reset the panel criterion back to Asian so the subsequent viewer steps start from a consistent
   baseline: in the Filter Panel, click the **Asian** row in the RACE filter card. The scatter-plot
   zoom is still applied, so the count does not come back to trueCount_panel here — wait for it to
   settle strictly BELOW trueCount_panel.

   Verify: with the scatter zoom still applied, the filtered row count is strictly BELOW
   trueCount_panel. Restoring the panel criterion undoes only the panel's own contribution — the
   viewer-driven narrowing is still in effect, and the count returns to trueCount_panel exactly only
   once the Scatter Plot is closed in Step 8.

7. Add a **Bar Chart** viewer. Right-click the Bar Chart and select **On Click > Filter**. Click a
   bar in the Bar Chart that corresponds to a non-empty subset of the currently filtered rows. Wait
   for the filtered row count to update.

   Verify: before the bar is clicked, the row count equals **trueCount_panel** exactly — the Scatter
   Plot has been closed, so the narrowing measured here is the bar chart's alone. Verify: the Bar
   Chart's **On Click** setting really reads Filter before the click is delivered.

   Verify: the number of rows passing the filter is strictly less than trueCount_panel (bar-chart
   click narrows the panel-filtered set further — AND intersection) and strictly greater than zero.
   Verify: the active-filter counter still reads **1** (the bar-chart click did not add a panel card).

8. Close the Scatter Plot viewer using its close button. Wait for the filtered row count to update.

   Verify: immediately BEFORE the close, the row count is strictly below trueCount_panel — the zoom
   narrowing that the close is supposed to release is proven to be in place, so "it came back" is
   falsifiable.
   Verify: the number of rows passing the filter equals **trueCount_panel** exactly (the scatter
   zoom narrowing is gone; the panel criterion alone survives). Verify: the active-filter counter
   still reads **1** and the RACE card still carries exactly the seeded category.

9. Close the Bar Chart viewer. Wait for the filtered row count to update.

   Verify: immediately BEFORE the close, the row count is strictly below trueCount_panel — the
   bar-chart narrowing is proven to be in place before the close is credited with releasing it.
   Verify: the number of rows passing the filter equals **trueCount_panel** exactly (the bar-chart
   narrowing is gone; the panel criterion is the only active filter). Verify: the active-filter
   counter still reads **1** and the RACE card still carries exactly the seeded category.

Expected:

- Step 4: entry count == trueCount_panel and the zoom-and-filter mode reads "filter by zoom"; then
  0 < filtered row count < trueCount_panel AND active-filter counter == 1 AND the RACE card still
  carries only the seeded category.
- Step 5: trueCount_after_change differs from the zoomed value in Step 3 AND trueCount_after_change
  is less than the Caucasian full-dataset count (zoom is still active) AND is above zero AND the
  active-filter counter == 1.
- Step 6: after restoring the seeded category the count is still < trueCount_panel (the zoom is
  untouched).
- Step 7: entry count == trueCount_panel and the On Click setting reads Filter; then
  0 < filtered row count < trueCount_panel AND active-filter counter == 1.
- Step 8: pre-close count < trueCount_panel; post-close filtered row count == trueCount_panel AND
  counter == 1 AND the RACE criterion unchanged.
- Step 9: pre-close count < trueCount_panel; post-close filtered row count == trueCount_panel AND
  counter == 1 AND the RACE criterion unchanged.

### Scenario 2: Layout round-trip with three filtering sources active (GROK-18281)

Steps:

1. Ensure at least one RACE categorical filter is active in the Filter Panel (from Setup), a Scatter
   Plot with filter-by-zoom is present with a zoom applied, and a Bar Chart with On Click > Filter
   is present with a bar clicked. Record the current filtered row count as **trueCount_combined**.
   Start RECORDING errors from this point on — every error balloon that appears and every uncaught
   page error, collected continuously until Step 5 reads them (see Automation notes). Do not take a
   count here: an error balloon fades on its own within a few seconds, so a count compared at the
   two ends of the window would miss an error raised inside it.
2. Save the current layout: open the **View** menu > **Layout** > **Save to Gallery**. Name the
   layout `probe-compose-viewer`.

   Verify: exactly ONE new layout of the current user appeared for this table. A save that produced
   none would leave Step 4 re-applying nothing, and a save that produced several would leave it
   re-applying an unknown one.
3. Perturb the layout: close the Filter Panel using its close button. Wait one second.
4. Re-apply the saved layout: open the layout gallery, locate `probe-compose-viewer`, and apply it.
   Wait 3 seconds for the filter cards to rebuild.
5. Verify: nothing was recorded on either error channel since Step 1 — no error balloon appeared at
   any moment and no page error fired. Verify: the Filter Panel closed in Step 3 is back with its
   cards rebuilt — at least one filter card is present inside the Filters viewer, read off the DOM
   and not through `getFiltersGroup()`, which creates a panel when none exists and would manufacture
   its own green — and the RACE criterion restored. Verify: the number of rows passing the filter is
   strictly below **rowCount** (the panel filter is active again; the exact count may differ from
   trueCount_combined because in-viewer filters restore asynchronously — assert only that the count
   is below rowCount, not an exact value).

Teardown:

- Delete the probe layout: in the layout gallery, select `probe-compose-viewer` and delete it
  (even if the scenario failed).
- Reset all filters: click the reset button in the Filter Panel header.

Expected:

- Step 5: nothing recorded on either error channel across the whole window; the Filter Panel is
  rebuilt with the RACE criterion restored; filtered row count < rowCount.

### Scenario 3: "Selection to filter" (Data menu command) — panel card state reflects surviving data (GROK-16713)

Steps:

1. Ensure the Filter Panel is open and a RACE categorical filter is active with only **Asian**
   selected. Record **trueCount_panel** (Asian rows only).

   Verify: no viewer other than the grid and the Filter Panel is present. A Scatter Plot or Bar
   Chart left over from an earlier scenario would add its own narrowing, and the row effect measured
   in Step 5 would no longer be readable — a zero row count could not be told apart from a dead
   command.

2. Select all currently visible rows: click in the grid and press **Ctrl+A**. Verify the grid's
   selection count matches trueCount_panel.
3. Select rows for a different single category (e.g. Caucasian) by Ctrl-clicking the **Caucasian**
   row in the RACE filter card. At this point the grid shows only Asian rows (the filter is still
   active) but the selection includes Caucasian rows — which are currently hidden by the filter
   (they are selected at the dataframe level but not visible in the filtered view).

   Verify: the second category really exists in the dataset (its full-dataset row count is above
   zero) and the selection grew to exactly trueCount_panel plus that count. The two categories are
   disjoint, so this pins both stages; without it the Step 5 arithmetic could be satisfied by zeros.

4. Invert the selection: open the top menu **Select** and choose **Invert**. Record
   **inverted_selection_count** — the number of rows selected after the inversion.

   Verify: inverted_selection_count equals the total row count minus the pre-inversion selection
   size exactly — the complement, not merely a different number — and lies strictly between zero and
   the total row count. Verify: the filter and the selection DISAGREE at this point (the filtered
   row count is trueCount_panel and differs from inverted_selection_count), so the Step 5 equality
   cannot be satisfied by a command that does nothing at all.
5. Run **Select > Selection to Filter** (top menu **Select**, leaf **Selection to Filter**). Wait
   for the filtered row count to update.

   Verify: the number of rows passing the filter equals **inverted_selection_count** (the filter
   was replaced by a verbatim copy of the selection — the selected rows are now the filtered rows).
   Verify: **Asian** — the panel card's own category — is now among the categories with zero
   surviving rows. That divergence between what the card is set to and what the data holds is the
   condition GROK-16713 was about, and it is what makes the card check below meaningful.
   Verify: the RACE card still has **Asian** as its criterion (the command copies the selection onto
   the filter; it must not rewrite the card's settings). Verify: the RACE card redraws — its canvas
   differs from the frame captured immediately before the command was issued, while that same canvas
   stayed pixel-identical across the idle window that preceded the command (GROK-16713 regression
   guard).

Teardown:

- Reset all filters: click the reset button in the Filter Panel header.
- Clear the selection: open the **Edit** menu and choose **Select None** (or press Escape while
  focus is on the grid).

Expected:

- Step 1: no viewer besides the grid and the Filter Panel is present.
- Step 3: the second category's full-dataset count is above zero and the selection equals
  trueCount_panel plus that count.
- Step 4: inverted_selection_count == rowCount − pre-inversion selection size, strictly between 0
  and rowCount, and different from the current filtered row count.
- Step 5: filtered row count equals inverted_selection_count, and Asian is now 100% filtered out.
  The RACE card keeps Asian as its criterion, and the card canvas redraws — different from the
  pre-command frame, after having been identical across the preceding idle window.

### Scenario 4: The remaining four viewer-driven filtering channels compose the same way

Scenario 1 covers two of the six ways a viewer can narrow the dataframe filter: the Scatter Plot's
zoom and the Bar Chart's click. This scenario covers the other four. They are kept apart from one
another rather than merged because each opens a different code path and, more importantly, each has
a different way of letting go again — which is the half that shows the narrowing really belonged to
that viewer.

Run the following four blocks in order. Each one starts from the same place: the RACE panel
criterion is the only thing filtering, so the row count equals trueCount_panel.

| Viewer | Filtering gesture | Revert gesture |
|---|---|---|
| Histogram | drag the max handle of the histogram's own range slider — painted on the viewer canvas, near its bottom edge — to the left across roughly 40% of the track | drag the same handle back to the right end of the track (its range saturates at the maximum there); double-clicking the slider does NOT restore the full range |
| PC Plot | drag the max handle of a per-axis slider downwards across roughly half its height | double-click that slider |
| Trellis Plot | set **On Click > Filter** from the viewer's context menu, then click in the LOWER part of a cell — the chart drawn in the middle of a cell swallows a click aimed at its centre | press Escape |
| Pie Chart | set **On Click > Filter** from the viewer's context menu, then click a segment | click outside the pie |

Steps, for each of the four:

1. Before adding the viewer, confirm the row count is trueCount_panel. If it is not, something from
   an earlier block is still filtering and nothing measured below would be attributable to this
   viewer — stop and report that rather than continuing.
2. Add the viewer and confirm it is present in the view.
3. Prepare it so its gesture reaches the filter. For the Histogram that means confirming it is set
   to filter by its range (a Histogram that only zooms would never move the row count, and the
   block would fail for the wrong reason). For the PC Plot it means turning its per-axis filter
   sliders on and confirming they really are on — with them off there is no slider to drag and the
   block would fail on the actuation rather than on the composition. For the Trellis Plot and the
   Pie Chart that means
   choosing **On Click > Filter** in the viewer's own context menu, and then confirming the
   viewer's On Click setting really reads Filter — driving the menu without checking would pass on
   a build where the menu item does nothing.
4. Perform the filtering gesture from the table above.

   Verify first that the gesture was DELIVERED, and stop with that as the failure if it was not: for
   the Histogram, that a range handle was found in the slider band at all (the canvas cursor reads
   `ew-resize` over it); for the Trellis Plot, that the click left a current cell behind. A gesture
   that was merely issued at a plausible place is not a delivered one, and the row-count checks below
   would otherwise report the miss as a filtering defect.

   Verify: the row count is strictly below trueCount_panel and strictly above zero — the viewer
   narrowed the already-filtered set rather than replacing it, and did not wipe every row. Verify
   the active-filter counter in the panel header still reads 1 — the viewer did not create a panel
   card. Verify the RACE card still holds exactly the category it was seeded with — the viewer did
   not rewrite the panel criterion.
5. Perform the revert gesture from the table above, leaving the viewer open.

   Verify: the row count is back to trueCount_panel exactly, and the RACE card still holds its
   category. Exactly this viewer's contribution was released and nothing else was.
6. Close the viewer.

   Verify: the row count is trueCount_panel and the counter reads 1.

Expected:

- Step 1 of each block: row count == trueCount_panel on entry.
- Step 4 of each block: the gesture is delivered (histogram handle found, trellis cell made current),
  then 0 < row count < trueCount_panel, counter == 1, RACE criterion unchanged.
- Step 5 of each block: row count == trueCount_panel with the viewer still open.
- Step 6 of each block: row count == trueCount_panel, counter == 1.

## Automation notes

Panel selectors, the card-repaint measures, the layout-apply readiness barrier and the console-error
rule are in `.claude/skills/grok-browser/references/viewers/filters.md`; the Histogram range slider
and its `filteringEnabled` arming are in `histogram.md`, the PC Plot axis handles and `showFilters`
in `pc_plot.md`, and trellis cell addressing with its delivery witness in `trellis_plot.md`, all in
the same `.claude/skills/grok-browser/references/viewers/` folder.

- **Setup — the RACE criterion is a DECLARED SETUP SUBSTITUTION**: it is applied through
  `fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']})`, and
  changed and reset the same way in Steps 5-6. The canvas row click is the equivalent UI path, but
  the panel criterion is never the gesture under test here — the viewer gestures are.
- **Both actuation paths were measured for the Step 5 criterion change, and they agree.** Step 5 is
  the one place where changing the criterion is itself the actuation rather than the scene, and a
  card click raises one notification `updateOrAdd` does not (`FILTER_CRITERIA_CHANGED` on the
  frame's event bus), so "the zoom survived" could in principle have been a property of the weaker
  path. Measured on dev 2026-08-18: with the scatter zoom already narrowing 354 rows to 222, a click
  on the RACE card's "Black" category-name row left 127 rows filtered where that category holds 157
  in full, with `zoomAndFilter` still reading `filter by zoom` — the zoom kept biting through a
  card-driven criterion change, which is what the step asserts through the API. The row order of the
  RACE card body is Asian, Black, Caucasian, Other, and a category-name click is an EXCLUSIVE
  select, so single-category criteria like these are reachable by gesture while the multi-category
  sets other scenarios need are not.
- **Scenario 1 Step 2**: a fresh Scatter Plot already comes up with `zoomAndFilter` reading
  `filter by zoom`, so the step asserts its OUTCOME by reading the property back before the drag
  rather than driving a menu to a value it already holds. A build that changes the default fails
  here first and by name.
- **Scenario 1 Step 7**: `onClick = Filter` and `splitColumnName = SEX` are set through the Bar
  Chart's properties and `onClick` is read back before the click — the gesture under test in this
  step is the BAR CLICK. Scenario 4 drives **On Click > Filter** through the real context menu for
  the Trellis and Pie charts, where the default is not Filter.
- **Scenario 1 step order**: the Scatter Plot close (Step 8) runs BEFORE the Bar Chart block
  (Step 7), so each viewer's narrowing is measured from the same panel-only base. The assertions
  are unchanged; only the run order differs from the numbering.
- **Scenario 1 Step 3**: the zoom is a real mouse drag (`page.mouse` move/down/move/up) on the
  Scatter Plot's `[name="canvas"]` / `[name="overlay"]`.
- **Scenario 1 Step 6**: the zoom is still applied here, so do NOT wait for trueCount_panel — wait
  for the count to settle strictly BELOW it.
- **Scenario 2 Step 1** spans ~25 s, longer than a balloon lives, so both error channels must
  accumulate across the whole window rather than be sampled at its ends; disconnect the observer in
  the step's `finally`.
- **Scenario 2 Step 3**: the panel must be closed by a real click on its title-bar close control. Do
  NOT substitute `dataFrame.resetFilter()` — that leaves the panel and all its cards standing, so
  the re-apply never rebuilds them and the GROK-18281 path is never entered.
- **Scenario 2 Step 5**: read the rebuild off the DOM (`[name="viewer-Filters"] .d4-filter`), never
  through `tv.getFiltersGroup()`, which creates a panel when none exists and would manufacture its
  own green.
- **Scenario 2 teardown**: delete the probe layout with `grok.dapi.layouts.delete(saved)` in a
  `finally` block, even on failure.
- **Scenario 3 Steps 2-3 — DECLARED SETUP SUBSTITUTION**: the two selection-building stages are
  placed through the dataframe — `df.selection.copyFrom(df.filter)`, then
  `df.selection.init(i => kept[i] || RACE[i] === 'Caucasian')` — because the UI equivalents need
  grid keyboard focus the headless run does not reliably hold. The gesture under test is the Step 5
  menu command, so both stages are setup, and each stage's resulting selection size is asserted
  before the command runs. This bullet is the sanction the spec cites.
- **Scenario 3 Step 4**: drive the REAL leaf with `driveTopMenuLeaf(page, ['Select', 'Invert'])`; no
  `selection.invert()` substitution — the menu path exists and the helper is already in the step.
- **Scenario 3 Step 5 — the card check is a REPAINT, not a count.** The gesture under test is the
  Step 5 menu command, and the only observable for "the card reflected the new filtering" is the
  card canvas repainting, so the assertion is honestly a floor and all three measures that make
  that floor discriminating are required here.
- **Scenario 4 viewer roots**: `[name="viewer-Histogram"]`, `[name="viewer-PC-Plot"]`,
  `[name="viewer-Trellis-plot"]`, `[name="viewer-Pie-chart"]`; each is closed through its own dock
  title-bar close control, like the Scatter Plot and Bar Chart in Scenario 1.
- **Scenario 4 — arm each channel and assert the property rather than assume it**, so a block that
  fails does so on the composition and not on a disarmed viewer. The Histogram's `valueColumnName`
  is set to AGE so the slider spans a real numeric range; Trellis and Pie start at `None` /
  `Select`, so their **On Click > Filter** is driven through the real context menu and the resulting
  `onClick` is read back.
- **Scenario 4 — Histogram**: the block asserts an EXACT return to trueCount_panel after the revert
  drag, so a drag that stops short fails instead of passing as near enough.
- **Scenario 4 — Trellis**: which cell to click is data-dependent — pick the one whose resulting
  count lands in the expected band rather than a fixed index, since several combinations are
  structurally empty and give 0, which is why the block asserts strictly above zero as well as
  strictly below trueCount_panel.
- **Readiness barriers**: wait for filter cards after opening the panel, and wait for the filtered
  row count to change after each viewer gesture before reading it.
