---
feature: trellisplot
realizes_atlas:
  - trellisplot.cp.split-and-pick-inner
  - trellisplot.int.split-columns-drive-inner-viewer-grid
  - trellisplot.int.viewer-type-change-control-panel-axes
realizes:
  - viewers.trellis-plot
priority: p0
target_layer: playwright
boot_lane: server
coverage_type: smoke
related_bugs:
  - id: github-964
    status: fixed
  - id: GROK-19633
    status: fixed
  - id: GROK-19890
    status: fixed
  - id: GROK-19634
    status: fixed
  - id: github-3015
    status: fixed
  - id: github-3067
    status: fixed
  - id: GROK-19673
    status: fixed
  - id: GROK-15494
    status: fixed
  - id: GROK-19902
    status: fixed
realized_as:
  - trellis-plot-split-and-pick-inner-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 1"
    expectation: >-
      The Trellis Plot viewer is present in the view exactly once immediately
      after it is added, and the dataset open plus the viewer add raise no
      uncaught page error and no non-benign console error (github-964 smoke
      guard). That floor is read on a window opened after login, so the
      platform's own bootstrap noise can neither satisfy nor fail it.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      The number of .d4-trellis-plot-cell DOM elements equals 8
      (xCategoriesCount=2 for SEX multiplied by yCategoriesCount=4 for RACE,
      both read from the TrellisPlotViewer JS object; the product equals the DOM
      cell count — no clamp fires because both axis counts are below 7.5).
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      Adding CONTROL as the second X-axis column through the (+) add-column
      control takes the .d4-trellis-plot-cell DOM count to 16
      (xCategoriesCount=4 for the SEX*CONTROL product, yCategoriesCount=4 for
      RACE) without any page reload or manual refresh, and the same actuation
      carries the GROK-19673 guard: hovering the CONTROL row previews (the grid
      moves off the committed 8 cells), Escape with the pointer still on the row
      reverts the grid to 8 cells and the X axis to its single committed column,
      and only the click commits both the second column and the 16-cell product.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      Removing the second X-axis column (CONTROL) through the blank first entry
      of that column's own selector list takes the X axis back to the single
      column SEX AND the .d4-trellis-plot-cell DOM count back to 8
      (xCategoriesCount=2, yCategoriesCount=4). Both witnesses are required: the
      cell count alone could return by coincidence, and only the column list
      says WHICH column went away. With Step 6 adding through the (+) control,
      the two steps close a full UI round-trip on the split-column selector —
      neither half is driven by a property write.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Each inner-type switch (Scatter plot, Bar chart, Histogram, Pie chart)
      fires the d4-trellis-plot-viewer-type-changed event with the type name,
      and after settle two independently sampled cells' canvases have
      non-identical pixel content — the GROK-19633 "same plot in every cell"
      failure class and the GROK-19890 blank-Sparklines class are both absent at
      the four exercised types.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      After the no-data excursion and the switch back to Pie chart (github-3015
      guard), each of the two sampled cells re-renders off the Bar chart content
      it held immediately before that switch rather than retaining it, and the
      two cells stay distinct from each other — the stale-render class does not
      reproduce, and two identically blank cells cannot pass as re-rendered.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      After hiding the control panel the type-selector control is no longer
      visible (zero-size under the hidden panel; the element may persist in the
      DOM), while a JS read of look.viewerType returns the same type that was
      active before hiding — the other-channel read confirms the setting was not
      lost.
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      After re-applying the saved layout the view's viewer set matches the saved
      set (both Trellis plots are back; the subsequently added Scatter plot is
      gone), and EACH Trellis plot returns with its OWN saved settings — one
      with xColumnNames ['SEX', 'CONTROL'], yColumnNames ['RACE'] and viewerType
      'Pie chart', the other with xColumnNames ['RACE'], yColumnNames ['SEX']
      and viewerType 'Bar chart'. Changing the inner viewer type of one restored
      Trellis plot afterwards leaves the other's inner type unchanged
      (GROK-15494 guard: Apply copies once and the two viewers stay
      independent). With the second Trellis plot closed again, the surviving one
      still reads the saved configuration and its .d4-trellis-plot-cell DOM
      count equals 16 (the restored 4*4 product under the viewport clamp).
  - anchor: "Scenario 3 Step 7"
    expectation: >-
      After setting On Click to Select and clicking the corner of one Trellis
      cell, the clicked cell announces its own column-to-category combination
      and df.selection.trueCount equals the (non-zero) number of dataset rows in
      that combination, and the viewer reads On Click = Select with Row Source =
      Selected; the subsequent ribbon Save completes with the project resolvable
      server-side by its typed name and raises no uncaught page error.
  - anchor: "Scenario 3 Step 11"
    expectation: >-
      After Close All and a real reopen of the saved project, the restored
      Trellis plot viewer is present in the view, xColumnNames equals ['SEX',
      'CONTROL'], yColumnNames equals ['RACE'], viewerType equals 'Pie chart',
      the restored category product reads 4 x 4 (xCategoriesCount 4,
      yCategoriesCount 4), props.rowSource equals 'Selected', props.onClick
      equals 'Select', and the reopen path raises no uncaught page error and no
      non-benign console error (GROK-19902 regression guard). The selection that
      was live when the project was saved does NOT come back:
      df.selection.trueCount reads 0 on the reopened view although Step 7 saved
      a non-zero selection [DOM 2026-08-12]. This step is the section's only
      reading of selection persistence, and it grades that loss explicitly
      against the count Step 7 recorded. Because Row Source = Selected plots the
      selected rows, the restored viewer consequently has nothing to plot and
      its .d4-trellis-plot-cell DOM count is 0; the empty grid is paired with a
      liveness witness — selecting every row fills it to the full 16 cells with
      the restored split columns unchanged. Any intermediate count in the
      liveness state means split columns were lost. The step then clears the
      selection again (the grid empties back to 0) and moves Row Source to All:
      with NOTHING selected the restored viewer must still paint the full 16
      cells with xColumnNames ['SEX', 'CONTROL'] — the proof that it assembles
      its grid without a selection at all, which the selection-driven witness
      cannot give.
  - anchor: "Scenario 3 Step 12"
    expectation: >-
      A second project saved with On Click = Select and Row Source = Selected
      and NO pre-selection (df.selection.trueCount is 0 at save time) reopens
      with the Trellis restored, raises no uncaught page error and no non-benign
      console error on the viewport-restore path, and the restored viewer reads
      Row Source = Selected with xColumnNames ['SEX', 'CONTROL'], yColumnNames
      ['RACE'] and the 4 x 4 category product. Its restored selection is empty
      because it was empty AT SAVE TIME — the state was carried faithfully, not
      lost the way Step 11's was — and that empty-at-save shape is the viewport
      the GROK-19902 crash lived on, which is what this step and only this step
      exercises. Its restored cell grid is EMPTY and that is the correct result
      — Row Source = Selected with nothing selected has no rows to plot — so the
      empty grid is paired with a liveness witness: selecting every row fills
      the grid to the full 16 cells while the restored split columns stay
      ['SEX', 'CONTROL'].
---

# Trellis Plot — Split columns, inner-type switching, and persistence (p0)

## Setup

1. Close all open views.
2. Open `System:DemoFiles/demog.csv` (the live DemoFiles demog golden dataset — its
   column cardinalities are: SEX=2, RACE=4, CONTROL=2 bool).
3. Add a Trellis Plot viewer to the table view (Add viewer > Trellis Plot, or from
   the toolbox icon).

## Scenarios

### Scenario 1: Split-column assignment drives cell count (github-964, GROK-19673)

Steps:
1. Verify the Trellis Plot viewer appears without any console error (github-964
   smoke guard — the viewer must be addable at all).
2. In the Trellis Plot control panel, open the X-axis column selector and select
   **SEX**.
3. Open the Y-axis column selector and select **RACE**.
4. Verify that the viewer grid shows 8 cells (2 X categories for SEX × 4 Y categories
   for RACE).
5. Click the **+** control next to the X-axis column selectors to add a second X
   category column. Hover over **CONTROL** in the list that opens — the grid may
   transiently preview the hovered column, but nothing is committed yet: press Escape
   with the pointer still on the row and the list closes, the grid falls back to 8 cells
   and the X axis reads **SEX** alone (GROK-19673 guard: hover previews, only click
   commits). Reopen the same **+** list and click **CONTROL** to commit it as the second
   X-axis column.
6. Verify that the viewer grid now shows 16 cells (4 X categories × 4 Y categories)
   without any page refresh or viewer recreation.
7. Open the column selector of the second X-axis column (**CONTROL**) and take the blank
   entry at the very top of the list that opens — that entry drops this one split column.
   Move the pointer straight onto the blank entry: the list previews the row under the
   pointer, so passing over the other rows on the way changes the split columns. Verify the
   X axis reads **SEX** alone again and the grid returns to 8 cells.

### Scenario 2: Inner-type switching renders per-cell data (GROK-19633, GROK-19890, GROK-19634, github-3015, github-3067)

Steps:
1. Confirm X = **SEX**, Y = **RACE** (2×4 = 8 cells from Scenario 1).
2. In the control panel, open the inner viewer type selector (the viewer-type
   dropdown at the top of the control panel).
3. Switch the inner viewer type to each of the four deep types in sequence — first
   **Scatter plot**, then **Bar chart**, then **Histogram**, then **Pie chart** —
   waiting for the cells to settle after each switch. After each switch, visually
   confirm that two independently chosen cells (e.g. top-left and row 2, col 2)
   display distinct, non-identical chart content.
4. From the Pie chart type, clear the Y-axis column by selecting the empty/none entry
   from the Y selector and switch the inner type to **Bar chart** — the viewer is now in
   a no-data configuration. Then put **RACE** back on the Y axis while still on Bar
   chart and let the cells settle: the Bar chart content they now hold is the reference
   the next step compares against. Finally switch the inner type back to **Pie chart**.
5. After the switch back to Pie chart in step 4, wait for the cells to settle and confirm
   that two independently chosen cells display the current type's content (Pie chart) and
   that each of them has moved off the Bar chart content it held immediately before the
   switch — no stale canvas from the previous type survives the excursion.
6. With the inner type set to **Pie chart**, set the inner viewer's **Marker Color**
   to **RACE** through the gear icon > inner viewer properties panel. Wait for settle
   and inspect the canvas of any non-empty cell. (Manual check — not separately
   automated: the color-coding bug classes ride the per-cell canvas checks of steps
   3 and 5.)
7. In the Trellis Plot control panel, toggle **Show Control Panel** off (or select
   the corresponding Context Panel option). Confirm the type-selector control is no
   longer visible in the viewer. Confirm the viewer still shows the type that was
   active before hiding.

### Scenario 3: Layout and project persistence (GROK-15494, GROK-19902)

Steps:
1. Set X = **SEX**, Y = **RACE**, add **CONTROL** as the second X column (4 X
   categories, 4 Y categories → 16 cells). Switch the inner type to **Pie chart**.
   Then add a **second** Trellis Plot viewer to the same view and give it deliberately
   different settings: X = **RACE**, Y = **SEX**, inner type **Bar chart**. Without a
   second viewer the layout round-trip cannot show whether the two keep independent
   settings.
2. Save the current view layout under a distinct name (e.g.
   "Trellis-split-test-layout") — see Automation notes for the channel this step is
   driven through and why.
3. Add a standalone **Scatter plot** viewer to the same view so the view now contains
   two viewers.
4. Re-apply the saved layout ("Trellis-split-test-layout") through the Layout menu.
5. Verify the grid cell count and the current X column assignment, Y column assignment,
   and inner viewer type as displayed in each Trellis Plot's control panel.
6. Verify that the view's viewer set matches the saved snapshot (the Scatter plot added
   in step 3 is gone) and that **both** Trellis Plots came back with their own saved
   settings — one with X = SEX + CONTROL, Y = RACE, inner type Pie chart; the other with
   X = RACE, Y = SEX, inner type Bar chart. Then switch the inner type of the first
   restored Trellis Plot to **Histogram** and confirm the second one still reads Bar
   chart (GROK-15494 guard: the two viewers must be completely independent); switch the
   first back to Pie chart. Close the second Trellis Plot and confirm the surviving one
   still shows X = SEX + CONTROL, Y = RACE, inner type Pie chart and 16 cells.
7. In the Trellis Plot property panel (gear icon > Context Panel), set **On Click** to
   **Select**, then click the corner of one Trellis Plot cell so the click registers a
   selection in the Table View, and note which column/category combination that cell
   stands for **and how many rows it selected** — step 11 compares against that number.
   Then, in the same property panel, set **Row Source** to **Selected**.
   Using the real ribbon **Save** button (File > Save), save the whole project to the
   server. Give it a distinct name (e.g. "Trellis-split-test-project").
8. Close all open views.
9. Reopen the saved project from the platform (browse to it in Files or Projects).
10. Wait for the Trellis Plot viewer to appear and settle.
11. Verify the current X column assignment, Y column assignment, inner viewer type,
    **On Click** and **Row Source** values as displayed in the Trellis Plot control panel,
    and confirm the browser console shows no uncaught exception from the reload path. Then
    count the rows selected in the reopened Table View: **none are selected** — the
    selection made in step 7 does not come back with the project [DOM 2026-08-12]. Compare
    that against the row count noted in step 7: this is the one place in the scenario where
    a selection is put through a save/reopen round-trip, so the comparison is the reading,
    not a side note. Because **Row Source = Selected** plots the selected rows, the restored
    viewer consequently draws **no cells**. Confirm it is idle rather than broken: select all
    rows in the Table View and watch the grid fill to 16 cells, still split by SEX + CONTROL
    over RACE. Then clear the selection — the grid empties again — and set **Row Source** to
    **All** in the property panel: with nothing selected at all the grid must come back to 16
    cells, still split by SEX + CONTROL over RACE. That is the separate question: the previous
    check shows the restored viewer draws when it is given a selection, this one shows it draws
    without one.
12. Repeat the persistence round-trip once more with **no pre-selection**: close all
    views, open demog.csv again, add a Trellis Plot with X = SEX + CONTROL, Y = RACE and
    inner type Pie chart, set **On Click** to **Select** and **Row Source** to
    **Selected** without clicking any cell first (nothing is selected in the Table View),
    save the project under a second distinct name (e.g. "Trellis-split-test-empty"),
    close all views, and reopen it. Verify the Trellis Plot is restored, its Row Source
    still reads Selected, the Table View still has nothing selected, and the console shows
    no uncaught exception from the viewport restore path. Here the empty selection is the
    state that was **saved** and came back intact — unlike step 11, where a real selection
    went in and nothing came back — and a project whose saved viewport is a Selected Trellis
    with no rows is the exact shape the reopen crash lived on. The restored viewer draws
    **no cells** — with Row Source = Selected and an empty selection there are no rows to
    plot. Confirm the viewer is nevertheless alive: select all rows in the Table View and
    watch the grid fill to 16 cells, still split by SEX + CONTROL over RACE.
13. After the scenario completes (whether it passed or failed), delete the saved layout
    "Trellis-split-test-layout" and both saved projects ("Trellis-split-test-project"
    and "Trellis-split-test-empty") from the server to avoid leaving orphaned state.

## Automation notes

- **CHANNEL — project save/reopen (Scenario 3 Steps 7-12) is UI-driven:** the persist goes through
  the ribbon Save flow (`saveProjectViaUI`), the API-assembled alternative not bringing the viewer
  back at all (refdoc: Layouts and projects). Reopen by id (`projects.find(id)` then `open()`),
  poll `grok.shell.tableViews` until the Trellis re-materializes, then read `xColumnNames` /
  `yColumnNames` / `viewerType` / `rowSource` / `onClick` off it.
- **CHANNEL — layout save/apply (Scenario 3 Steps 2-6) is the API on purpose:** `tv.saveLayout()`
  + `grok.dapi.layouts.save()`, re-applied with `layouts.find()` + `tv.loadLayout()`. The
  grok-browser skill ("Saving and restoring a layout") prescribes `grok.dapi.layouts` over the SAVE
  ribbon / Ctrl+S, which do not guarantee the saved layout can be found again, and no ribbon
  Save-Layout affordance for this flow is documented. The PROJECT half is the opposite case: there
  the ribbon Save IS the contract under test.
- **CHANNEL — Scenario 2 Step 4** runs on the props channel; the UI-driven switch with its event
  assertion is Step 3's subject.
- **WITNESS — cell count:** `.d4-trellis-plot-cell` under the viewer root, never document-wide,
  cross-referenced with `xCategoriesCount` / `yCategoriesCount` off the viewer handle.
- **WITNESS — event capture:** `viewer.onEvent('d4-trellis-plot-viewer-type-changed')` subscribed
  BEFORE each switch, the switch driven through the control-panel selector — the only channel that
  raises it (refdoc: Events). Assert it fired; the argument is the type-name string.
- **WITNESS — per-cell canvas diff:** `.d4-trellis-plot-cell canvas` → `toDataURL()`, compared by
  hash, taken AFTER a settle gate (100 ms render + 300 ms refresh debounce; poll the type-changed
  event or wait ≥ 500 ms). Never a global non-white pixel count.
- **WITNESS — canvas-only inner cells:** inner-cell content has no DOM signal beyond the canvas, so
  the assertion is a canvas delta (two cells differ), never absolute pixel values.
- **WITNESS — GROK-19673 (hover vs commit):** grade COMMIT SEMANTICS, not render stillness —
  Escape must revert the grid to the committed product, a click must commit the new column and the
  cell count (refdoc: The (+) picker's commit semantics). The transient rebuild's ABSENCE is a soft
  preview-regression signal asserted separately. The guard runs on the ADD path (refdoc: Add X
  column), so it and the Scenario 1 16-cell round-trip share ONE actuation.
- **WITNESS — split-column removal (Scenario 1 Step 7):** two UI routes with different granularity
  (refdoc: Removing a Split Column). Step 7 drops ONE column, so it takes the blank picker row;
  `Reset X columns` is the known alternative for a future whole-axis-reset scenario and is NOT
  interchangeable here.
- **WITNESS — the bound picker is never scanned** (refdoc: Removing a Split Column, pitfall 21):
  the pointer goes STRAIGHT to the blank row in one `page.mouse.move`. Step 7 therefore cannot
  reuse the tooltip scan Step 6 runs in the (+) picker, whose rows commit nothing — the two pickers
  look alike and the rule is asymmetric.
- **WITNESS — locating the blank row (Scenario 1 Step 7):** the row has no DOM node, so it is
  addressed through the popup's layout model (refdoc: Locating a row in the canvas picker). The
  spec re-checks the height identity BEFORE the pointer moves, failing with the measured numbers,
  and re-reads the split columns after the hover but BEFORE the click, so a mis-targeted row is
  dismissed with Escape rather than committed. The CONTROL slot is located structurally (refdoc:
  Locate the split slot STRUCTURALLY), never by position.
- **WITNESS — where the Scenario 2 Steps 4-5 measured window sits:** the graded transition is
  Step 4's LAST switch, Bar chart → Pie chart, with the Y column already restored. The github-3015
  class needs a per-cell baseline taken while the cells hold the intermediate Bar chart frame,
  which requires data on the axes. The no-data excursion (`yColumnNames = []` on Bar chart) is a
  state the viewer must SURVIVE and sits ABOVE the measured window — a hash there would have no
  canvas to read. It visits no Scatter plot: Step 3 already canvas-checks Scatter.
- **WITNESS — the second Trellis plot in the layout probe (GROK-15494):** two trellises with
  deliberately different settings, both under the clamp — a single-viewer probe cannot see the bug,
  one shared look satisfying it just as well. Grade order-independently, sorting the pair by inner
  type (refdoc: Layouts and projects — restore order is not guaranteed), and prove independence by
  MUTATING one restored viewer and reading the OTHER. Close the second before the project half: the
  reopen path reads "the" restored trellis.
- **WITNESS — reopen console/pageerror floor:** the build's two benign errors on this path are
  filtered (refdoc: Layouts and projects); the GROK-19902 guard is the UNCAUGHT pageerror floor on
  the reopen path, not the console floor on the save step.
- **WITNESS — the reopen cell grid is graded by DATA STATE, never by box size** — the empty grid it
  can land on is correct by design (refdoc: `Row Source = Selected`: the empty grid and the lost
  selection). Grade the DOM count against `rowSource` plus the selection count sampled at the SAME
  moment: `Selected` + empty owes an empty grid, every other state owes the full 4 × 4 = 16, and an
  intermediate count means split columns were lost. `xCategoriesCount` / `yCategoriesCount` read
  4 × 4 in both states and stay the config witness.
- **WITNESS — the selection does NOT survive the project round-trip**, by design (refdoc:
  `Row Source = Selected`: the empty grid and the lost selection). Step 11 asserts the LOSS against
  the save-time count carried over from Step 7, never `0 == 0`.
- **WITNESS — what separates Step 11 from Step 12** is the transition covered and the probe run
  afterwards, never the restored state, identical in both. Step 11 saves a LIVE selection and reads
  its loss, then clears the selection, switches `rowSource` to `'All'`, requires 16 cells with
  `xColumnNames` unchanged, and re-reads `onClick`. Step 12 saves an EMPTY selection — the viewport
  shape GROK-19902 crashed on — and stays on the selection channel.
- **WITNESS — Scenario 2 / Step 7 control-panel visibility:** after Show Control Panel off, the
  type selector is graded by computed zero size, never by DOM absence (refdoc: Hiding the control
  panel keeps the node), with `look.viewerType` as the other-channel witness. Scenario 3 Steps 5
  and 11 read `xColumnNames` / `yColumnNames` the same way.
- **Viewport-clamp caveat:** cell-count equality holds only below the per-axis clamp (refdoc: The
  category viewport clamp), which the prescribed SEX=2, CONTROL=2, RACE=4 stay under. At the clamp,
  equality belongs to cp.scroll-categories.
- **Empty-combination caveat (Scenario 1):** the structurally empty SEX/CONTROL/RACE combinations
  do NOT move the 16-cell count (refdoc: What packing actually collapses, pitfall 19). For
  Scenario 2's per-cell diffs, sample populated CONTROL=false columns to guarantee non-empty
  canvases.
- **Persistence cleanup:** the layout delete and BOTH project deletes run in `finally`, so a
  mid-scenario failure cannot orphan server state.
- **Spec must keep** (guard-rails against re-authoring regressions):
  - Scenario 1 Step 1's error floor is read IMMEDIATELY after the viewer is created, not only at
    the Step 7 / 11 / 12 reopen floors — else an error raised while ADDING the viewer is never
    graded and the github-964 guard stays half unasserted against the step's own promise.
  - That read is WINDOWED (both collectors sliced at marks taken after login and before the build,
    console keeping `isBenignError`), so it fails on the add rather than on bootstrap noise.
  - Step 1 also asserts the viewer root is present exactly once — the "addable at all" half.
  - Scenario 1 Step 6 stays UI-driven on the (+) `[name="add-x-column"]`, with the GROK-19673 guard
    on THAT SAME actuation — never split across channels, never substituted by an `xColumnNames`
    write.
  - It opens with a real mouse down/up and locates CONTROL by real-mouse scanning the UNFILTERED
    list through the popup's column tooltip.
  - NEVER type into the search box: filtering is not part of the hover-preview contract and
    corrupts the Escape semantics under test.
  - Because the preview writes through (refdoc: The (+) picker's commit semantics), the dwell
    asserts only that it is ALIVE (grid off the committed 8 cells, or a second column bound);
    NON-COMMIT is proven by Escape with the pointer still on the row — popup closes, grid back to
    8, `xColumnNames` back to `['SEX']`. The real-mouse click then asserts the commit
    (`['SEX', 'CONTROL']`, 4 × 4, 16 cells).
  - The `finally` dismisses a surviving popup but must NOT write the split columns back: that
    write is the actuation under test.
  - Scenario 1 Steps 6 and 7 are ONE closed UI round-trip on the split-column selector — add
    through the (+), remove through the BLANK first row of that column's own picker. Neither half
    may become a `props.xColumnNames` write: this pair is the section's only two-directional
    exercise of the split-column UI.
  - Step 7 grades on TWO mandatory witnesses, the list back to `['SEX']` AND the count back to 8 —
    a returning cell count alone cannot say WHICH column went away.
  - The picker opened from an ALREADY-BOUND selector is NEVER scanned — the preview destroys the
    committed slot on every row it crosses (refdoc: Removing a Split Column, pitfall 21), so a
    scan silently empties the axis.
  - Step 7 moves the pointer ONCE onto the blank row derived from the layout model, checks the
    removal BEFORE clicking, and on a miss dismisses with Escape and fails with the measured
    geometry — no second candidate, no positional nudging, no arithmetic fallback.
  - `Reset X columns` must NOT be substituted into Step 7: it clears the WHOLE axis, a different
    product operation from dropping one column.
  - Scenario 2 Step 3 asserts the `d4-trellis-plot-viewer-type-changed` argument for every
    inner-type switch, driven through the control-panel selector — the event does not fire on the
    props channel, and a props set-then-read echo does not prove the switch was processed.
  - Canvas checks stay per-cell differentials (two populated cells differ; Step 5 compares against
    the intermediate Bar-chart frame), never a global non-white count.
  - The Step 4 baseline is taken AFTER the Y column is back and while the type is still Bar chart
    — moving the measurement into the no-data excursion, which has no canvas to hash, retires the
    guard silently.
  - Scenario 2 Step 7 asserts the hidden selector by computed zero-size visibility, never by DOM
    absence (the node survives the hide), plus the `look.viewerType` other-channel read.
  - It asserts the selector VISIBLE with the panel still on FIRST — a renamed or missing node reads
    zero-size in BOTH states, and the hide then passes on a control that was never there.
  - Scenario 3 Step 6 keeps the SECOND Trellis plot in the layout probe — the GROK-15494 guard is
    unconditional, and a single-viewer round-trip cannot fail on the bug. Both viewers' restored
    settings are compared against save time, order-independently; independence is proven by
    mutating one restored viewer and reading the OTHER; the second is closed before the project
    half.
  - Scenario 3 Step 6 persists the layout through `grok.dapi.layouts`, so the scenario text must
    not promise a ribbon Save-Layout button. The ribbon is required for the PROJECT half only.
  - Scenario 3 Step 7 grades the corner click on the EXACT row count of the clicked combination:
    capture `args.matchCondition`, count the matching dataset rows, assert that count non-zero and
    `selection.trueCount` equal to it — a bare "greater than 0" passes on any cell and hides an
    off-target click.
  - Scenario 3 Steps 7-11 keep the FULL project round-trip: real ribbon Save, Close All, reopen by
    id, then reads of `xColumnNames` (`['SEX', 'CONTROL']`), `yColumnNames` (`['RACE']`),
    `viewerType` (`'Pie chart'`), the 4 × 4 cell-grid product graded by the restored DATA STATE
    (never an intermediate count, never an unconditional `[0, 16]`), the selection count read at
    the same moment and graded against Step 7's, `onClick` (`'Select'`), `rowSource`
    (`'Selected'`), plus the zero-uncaught-page-error floor on the reopen path.
  - A "Save dialog opened and no new warning appeared" floor must never be substituted for that
    round-trip: the deterministic reopen-and-read is driveable here and the floor version does not
    guard the regression.
  - Scenario 3 Step 12 keeps the no-pre-selection variant (empty selection with Row Source =
    Selected) — that is the viewport-restore shape GROK-19902 crashed on, which the pre-selected
    variant does not cover.
  - Step 12 rebuilds a project-less workspace first so the ribbon Save creates a SECOND project
    rather than re-saving Step 11's, and re-asserts the reopened project still carries an EMPTY
    selection — a reappeared selection would silently retire the shape.
  - An EMPTY cell grid under Row Source = Selected with an empty selection is CORRECT and is
    asserted as such — not tolerated, not hedged (confirmed by design; no ticket exists and none
    will be filed).
  - Two rewrites are forbidden there: widening back to an unconditional `[0, 16]`, which admits the
    very regression the step catches, and re-tying the expectation to the viewer's box size in any
    form, measured false (refdoc: Layouts and projects records that reading as withdrawn).
  - The grading input is `rowSource` + selection count sampled together with the DOM count, and the
    empty branch ALWAYS carries the liveness witness (`df.selection.setAll(true)` must fill the
    grid to 16 with `xColumnNames` unchanged) — else "correctly empty" passes as "restored dead".
  - The 0-cell read comes AFTER the settle loop polling for 16 cells has run OUT, never straight
    after the viewer poll — an unpainted grid reads 0 too, so a grid that does paint passes green.
  - On the reopen path the pair is a constant, not a fork: both saved projects come back `Selected`
    with an empty selection, so the shared helper carries the empty case only and re-checks the
    state instead of branching. A "restored selection is non-empty" branch would be unreachable and
    must not return as a hedge.
  - "The restored selection paints all 16 cells" must never be reinstated in any form: it describes
    a state this pair does not reach.
  - A non-zero restored selection count must FAIL rather than be adapted to — it means the product
    behaviour CHANGED and both steps need re-deriving.
  - Steps 11 and 12 are told apart by the TRANSITION each covers and the PROBE each runs, never by
    the restored state, identical in both. Neither may be dropped as a duplicate, and neither may
    be given a fabricated difference in the restored state to look distinct.
  - Scenario 3 Step 11 ends with the ROW-SOURCE probe: clear the selection (the grid must empty
    back to 0 — the anti-vacuity half, the preceding witness having left every row selected), set
    `rowSource = 'All'`, require exactly 16 cells with `xColumnNames` unchanged. Not a second
    liveness witness: that one proves the viewer draws WITH a selection, this proves it draws
    WITHOUT one.
  - That probe stays AFTER the `savedSelectionCount` / `r.selected === 0` comparison, or it masks
    the reading the step exists for, and is deliberately NOT replicated into Step 12, where the
    saved `Selected` wiring is the subject.
  - `reopenAndReadTrellis` prints one non-assertive line with `rowSource`, the selection count and
    the cell count after each reopen — else the numbers that decide the grading appear in no
    assertion message on a passing run.
  - Restored On Click = Select is not re-exercised after the reopen; its effect is proven pre-save
    at Step 7.
  - The Scenario 3 probe deletes (saved layout and BOTH saved projects) run in a `finally`, so a
    mid-step failure cannot orphan them on the server.
