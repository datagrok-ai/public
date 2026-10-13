---
feature: trellisplot
realizes_atlas:
  - trellisplot.cp.scroll-categories
  - trellisplot.int.pack-categories-vs-viewport-scroll
realizes:
  - viewers.trellis-plot
priority: p2
target_layer: playwright
boot_lane: local
coverage_type: regression
realized_as:
  - trellis-plot-scroll-categories-spec.ts
related_bugs: []
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      After clicking the '+' X-axis icon once, the .d4-trellis-plot-cell DOM
      count increases by exactly one row of cells (equal to yCategoriesCount),
      confirming one additional X-axis category entered the viewport.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After clicking the '-' X-axis icon once, the .d4-trellis-plot-cell DOM
      count returns to its pre-expansion value, confirming the DOM-count
      round-trip is symmetric.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      When the viewport shows the minimum number of X-axis categories (one), the
      '-' axis icon is disabled, and a real click on it neither pages the grid
      nor is accepted at all, while the '+' icon beside it still responds; when
      all X-axis categories are visible, the same holds with the roles swapped.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      No uncaught page error and no console error is raised by the paging in
      Steps 4-6, read as a delta over that window rather than as a global
      empty-collector assertion.
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      Clicking the X column selector's context-menu "Reset X columns" item while
      categories are paged clears the X assignment (leaving one category row of
      cells) with no uncaught JavaScript error in the browser console.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      With Pack Categories ON, clicking the '+' X-axis icon increases the
      .d4-trellis-plot-cell count by exactly one row of cells (the number of Y
      categories); and once the viewport is opened all the way, the packed grid
      holds at most one cell row per POPULATED RACE/SEVERITY combination, which
      is strictly fewer cells than the full combination product would occupy.
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      Toggling Pack Categories OFF and opening the viewport all the way again
      restores the empty-combination cells: the cell count equals the full RACE
      x SEVERITY combination product times the number of Y categories, and is
      strictly greater than the packed count measured at the same fully opened
      viewport in Step 4. That difference is the packing effect, and the toggle
      itself raises no uncaught console error.
---

# Trellis Plot — Category viewport paging (+/- icons and pack-categories coupling)

## Setup

1. Close all open views.
2. Open `System:DemoFiles/demog.csv` (the live DemoFiles demog golden dataset —
   cardinalities relevant here: DIS_POP = 6, RACE = 4, SEVERITY = 5, SEX = 2).
3. Add a Trellis Plot viewer to the table view (Add viewer > Trellis Plot, or from
   the toolbox icon).
4. In the Trellis Plot control panel, set the X-axis column to **DIS_POP** and the
   Y-axis column to **RACE**.

## Scenarios

### Scenario 1: + / - axis icons page categories in and out (DOM-count round-trip)

At entry every DIS_POP category already fits the viewport, so the '+' icon starts
disabled and the '-' icon enabled — room for '+' has to be made by clicking '-'
first. The '+' and '-' icons on the X axis add or remove one category column from
the viewport.

Steps:
1. Confirm the Trellis Plot control panel shows X = **DIS_POP** and Y = **RACE**.
2. Count the number of Trellis Plot cells currently rendered in the grid (baseline
   cell count = visible X categories × the number of RACE categories).
3. Confirm the '+' X-axis icon is styled as disabled at entry and the '-' icon is
   enabled. Click the '-' X-axis icon once to page one category out; the '+' icon
   becomes enabled.
4. Click the '+' X-axis icon once and wait for the grid to settle. Count the
   Trellis Plot cells again.
5. Click the '-' X-axis icon once and wait for the grid to settle. Count the
   Trellis Plot cells again.
6. Shrink the viewport to its minimum by clicking '-' until the icon becomes
   disabled. There, click the disabled '-' icon as a user would and confirm nothing
   happens, then click the '+' icon beside it and confirm it still pages. Expand the
   viewport to its maximum by clicking '+' until that icon becomes disabled, and
   repeat the same pair of clicks with the roles swapped.
7. Confirm no uncaught JavaScript error appeared in the browser console during
   steps 3 through 6.
8. With the viewport in a mid-range paged state (not all categories visible), open
   the X-axis column selector, right-click the selector list header or use the
   context-menu Reset item to reset the X column assignment, and confirm the
   operation completes without a console error.

### Scenario 2: Pack Categories couples with viewport scrolling

Packing drops whole X-axis category combinations that no row carries, so the packed
grid is narrower than the full combination product; the +/- paging range adjusts
accordingly. The pair of columns matters: the collapse is only observable on an
X axis that actually has empty combinations. **RACE** combined with **SEVERITY**
has them on demog (Asian + Critical among others), whereas SEX combined with
SEVERITY has every combination populated and would collapse nothing.

Steps:
1. Set the X axis to two columns — **RACE** first, **SEVERITY** second — and the Y
   axis to **SEX**.
2. Note the current cell count in the grid (baseline) and confirm the viewer is
   paged (the '+' X-axis icon is enabled, so not every X-combination is on screen).
3. In the Trellis Plot control panel, toggle **Pack Categories** ON (if not
   already on). Wait for the grid to settle.
4. Click the '+' X-axis icon once and count the Trellis Plot cells. Then keep
   clicking '+' until it becomes disabled — the whole X axis is now on screen —
   and count the cells again.
5. Confirm the fully opened packed grid holds no more than one cell row per
   populated RACE/SEVERITY combination, and fewer cells than the full combination
   product would occupy.
6. Toggle **Pack Categories** OFF, wait for the grid to settle, open the viewport
   all the way again with '+', and count the cells.

## Automation notes

- **CHANNEL — two click channels on purpose** (refdoc: Category paging icons): the paging loops
  use `page.evaluate(el.click())`, honest only on icons already asserted enabled; Scenario 1
  Step 6's extremes use a REAL Playwright click, whose refusal IS the inertness claim, kept as a
  differential against the enabled sibling in the same state.
- **CHANNEL — props for preconditions:** `packCategories` and every axis-column assignment
  (Setup step 4, Scenario 1 step 8 setup and restore, Scenario 2 step 1) are prop writes, being
  settings rather than the tested action; the control-panel checkbox and column selector are the
  manual path. Exception: the RESET in Scenario 1 step 8 IS the tested action and goes through
  the real context menu.
- **CHANNEL — no-error floors:** `pageerror` + `console(error)` subscribed before login, benign
  classes filtered, each floor a DELTA over its own window (Step 7 over Steps 4-6, Step 8 across
  the reset, Scenario 2 Step 6 across the toggle).
- **WITNESS — cell count:** `.d4-trellis-plot-cell` scoped to the viewer root, never
  document-wide. The rendered count is `viewport.width × viewport.height`, so it equals
  (visible X) × (visible Y) only at a fully opened viewport — hence every packing comparison is
  taken after '+' has run to the disabled state.
- **WITNESS — +/- icons** (refdoc: X axis +/- icons, Category paging icons): located by their
  `name=` inside the axis container, never by child position; their disabled state is read as
  computed `pointer-events`, never `opacity`.
- **WITNESS — packing** (refdoc: What packing actually collapses): the distinct-(RACE,
  SEVERITY)-pair count is an UPPER bound, so the packed grid is graded `<=` and exact equality is
  taken only unpacked. The empty-combination precondition is computed from the live frame and
  asserted, so a demog that ever fills every combination fails loudly instead of turning the step
  into a no-op.
- **SCOPE:** this pair stays in the paging regime, well under the display limit (refdoc:
  pitfall 5) — the fully opened unpacked grid is 40 cells; the 250-cell and 100 M-combination
  guards belong to the atlas `edge_cases[]`.
- **Spec must keep (checklist for any future re-author of this pair):**
  - Scenario 2 stays on RACE + SEVERITY (refdoc: pitfall 19) — SEX × SEVERITY is fully populated
    and could not tell a working `packCategories` from a dead one.
  - Scenario 2 asserts the empty-combination precondition BEFORE the counts — else packing
    collapses nothing and the counts prove nothing.
  - Scenario 2 grades packed vs unpacked at the SAME fully opened viewport — else the two counts
    are read through different windows.
  - Scenario 2 Step 6 asserts the carried `packedMaxCells` > 0 BEFORE comparing it — `softStep` can
    swallow Step 4, and the sentinel -1 passes `toBeGreaterThan`, comparing the delta with nothing.
  - Scenario 2 never uses `xCategoriesCount` as the packing signal: it is packing-blind (refdoc:
    pitfall 26), so the assert could not fail.
  - No global `expect(pageErrors).toEqual([])`; every floor stays a windowed delta through the
    benign filter — else ambient bootstrap noise makes it false-red.
  - Scenario 1 Step 6 grades both extremes with a REAL Playwright click, as a differential
    against the enabled sibling; a synthetic click fires through `pointer-events: none` and
    proves no inertness. The `page.evaluate` click stays for the enabled-icon loops only.
  - The prose keeps the observed entry state (all DIS_POP categories fit, '+' starts disabled,
    '-' has to be clicked first) AND the spec asserts it before the first click — else the claim
    is assumed rather than graded.
  - Scenario 1 Step 8 keeps X on TWO columns (SEX + DIS_POP) and asserts the paged state before the
    reset — with the single Setup column nothing pages and the crash-prone path is never entered.
  - Scenario 1 Step 8 reads the menu only inside `.d4-menu-popup` and gates its synthetic item click
    on a non-zero box measured FIRST (refdoc: pitfall 22) — unscoped it grades the wrong menu, and
    ungated a synthetic click fires through `display: none` onto an item no user can reach.
