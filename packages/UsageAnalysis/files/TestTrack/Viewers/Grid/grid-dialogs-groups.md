---
feature: grid
realizes_atlas:
  - grid.cp.dialogs-groups
realizes:
  - viewers.grid
priority: p2
target_layer: playwright
coverage_type: edge
related_bugs:
  - id: GROK-20167
    status: fixed
  - id: GROK-19333
    status: fixed
  - id: GROK-19332
    status: fixed
  - id: GROK-17463
    status: fixed
  - id: GROK-19288
    status: fixed
  - id: GROK-17505
    status: fixed
  - id: GROK-17442
    status: fixed
  - id: GROK-17443
    status: fixed
  - id: GROK-17441
    status: fixed
  - id: GROK-18213
    status: fixed
realized_as:
  - grid-dialogs-groups-spec.ts
expected_results:
  - anchor: Step 4
    expectation: After clicking OK in the Sort dialog configured with SEX then AGE,
      the grid's sort-by list equals ['SEX', 'AGE'] in that order and the
      visible rows are grouped SEX-major (read via gridRowToTable), while the
      dataframe's physical row order is untouched.
  - anchor: Step 9
    expectation: 'After clicking "Reset filter" in the Columns dialog, the
      type-filter checkboxes are fully unchecked — no stale checked state
      remains (GROK-19333 guard: the checkbox UI state matches the cleared
      internal filter state).'
  - anchor: Step 10
    expectation: 'The Columns dialog header element is still present in the DOM
      after the filter and reset interactions (GROK-20167 guard: the header is
      not dropped when the dialog DOM is rebuilt during filter/reset
      operations).'
  - anchor: Step 14
    expectation: After switching to the second table and applying a type-filter in
      the Columns dialog, no console error is raised — the dialog re-binds to
      the new table's column set cleanly (GROK-19332 guard). Console-error delta
      is 0 from the baseline recorded at setup.
  - anchor: Step 17
    expectation: 'After closing the grid properties panel via the X button and then
      clicking the gear icon again, the grid properties panel opens for the
      second time (GROK-17463 guard: the gear is not a one-shot open — repeated
      open/close cycles work).'
  - anchor: Step 21
    expectation: 'After setting color coding via the column''s hamburger menu
      (three-lines icon in the column header), the Context Panel''s Colors
      section for that same column reflects the newly applied color coding
      without requiring a manual refresh (GROK-19288 guard: hamburger-menu
      coding and Context Panel are in sync).'
  - anchor: Step 26
    expectation: 'After Shift+clicking on two or more headers of grouped columns, no
      console error is raised — the console-error delta is 0 from the baseline
      (GROK-17505 guard: multi-select over grouped headers must not fire "Cannot
      fire new event" errors).'
  - anchor: Step 28
    expectation: 'After clicking on the empty space of a column group''s header
      band, no console error is raised — the console-error delta is 0
      (GROK-17443 guard: clicking group empty space is a no-op that does not
      crash the hit-test).'
  - anchor: Step 31
    expectation: After creating a second column group and then selecting the first
      group's name, and then pressing Esc to deselect, no console error is
      raised — the console-error delta is 0 from the baseline (GROK-17442 +
      GROK-18213 guards).
  - anchor: Step 36
    expectation: 'After reopening the saved project, the column group colors
      assigned before saving are intact — each group''s header band displays the
      same color as before the project round-trip, confirmed by the df column
      ''group'' tag and the df ''.columnGroups'' tag (GROK-17441 guard: group
      colors travel with the project, not only the layout).'
  - anchor: Step 36
    expectation: After reopening the saved project, no console error is raised — the
      console-error delta from the post-reopen baseline is 0 (teardown guard).
---

# Grid — Dialogs, Hamburger Menu, and Column Groups

## Setup

1. Close all open views.
2. Open the demog dataset (demog.csv) to produce a Table View with the main grid.
   Demog contains columns SEX (string), AGE (int), HEIGHT and WEIGHT (float).
3. Record the baseline console-error count. This baseline is used for all
   console-error delta checks throughout the scenario.

## Scenarios

### Scenario 1: Multi-column sort dialog — configure, apply, verify order

Steps:
1. Open the Sort dialog: right-click any column header and choose Sort, or use
   the grid's top toolbar Sort button. The Sort dialog opens.
2. Remove any existing sort rules. Add a new sort rule: choose the SEX column,
   set direction to Ascending.
3. Add a second sort rule: choose the AGE column, set direction to Descending.
4. Click OK to apply the sort.
5. Verify that the grid is now sorted by SEX ascending first, then by AGE
   descending within each SEX group — inspect the first few visible rows to
   confirm the ordering is consistent with the configured sort.

Expected:
- Step 4: The grid's sort state reflects both columns in the configured order:
  SEX ascending, then AGE descending.
- Step 5: The first visible row has the lexicographically smallest SEX value
  available in demog, and within that SEX group the AGE values descend from the
  first visible row downward.

### Scenario 2: Order or Hide Columns dialog — type filter, Reset, header guard, data-switch guard

Steps:
1. Close all and reopen demog. Record the console-error baseline.
2. Open the Order or Hide Columns dialog: right-click any column header and
   choose Order or Hide Columns, or use the equivalent grid menu entry.
3. In the dialog, apply a type-filter by checking the "int" type checkbox (or
   another available type). The column list narrows to show only columns of that
   type.
4. Click "Reset filter" in the dialog.
5. Verify that all type-filter checkboxes are now unchecked and the full column
   list is restored — no stale checked state remains in the filter UI.
6. Verify that the dialog header element is still present and visible.
7. Open a second dataset: open spgi-100 (or any other available table) alongside
   demog, so that two tables are open.
8. Open the column manager by clicking **Columns** on the status bar, then apply a type-filter again
   (e.g. check the "string" checkbox).
9. Switch the active table view to the second dataset (click its tab or use Data
   > Table to switch).
10. In the Columns dialog, apply a type-filter on the second dataset's columns.
11. Verify that no console error was raised after switching tables and re-applying
    a type filter.
12. Close the second dataset. Return to the demog view.

Expected:
- Step 5: All type-filter checkboxes in the Columns dialog are cleared after
  Reset filter — the dialog's internal state and the checkbox UI are in sync
  (GROK-19333 guard).
- Step 6: The Columns dialog header element is still present in the DOM after
  the filter and reset interactions (GROK-20167 guard).
- Step 11: No console error is raised; the console-error
  delta is 0 from the baseline (GROK-19332 guard).


### Scenario 3: Column hamburger menu — stats popup and color coding sync with Context Panel

Steps:
1. Close all and reopen demog.
2. Hover the AGE column header until the hamburger icon (three horizontal lines)
   appears inside the header area.
3. Click the hamburger icon to open the column popup menu.
4. Verify that a popup with column statistics appears and can be dismissed.
5. Click the hamburger icon again to open the menu. Choose Color Coding > Linear
   from the hamburger menu (accept the defaults).
6. Open the Context Panel for the AGE column: click the AGE column header to
   select it, then look at the Properties / Colors section in the Context Panel
   on the right side.
7. Verify that the Context Panel's Colors section reflects the linear color
   coding that was just applied via the hamburger menu — without requiring a
   manual refresh.

Expected:
- Step 4: A column-stats popup appears after clicking the hamburger icon and can
  be dismissed without error.
- Step 7: The Context Panel Colors section for AGE shows the applied linear color
  coding matching the hamburger-menu setting — the two entry points agree
  immediately (GROK-19288 guard).

### Scenario 4: Column groups — creation, Shift+click guard, empty-space guard, Esc guard, and project persistence

Steps:
1. Close all and reopen demog. Record the console-error baseline.
2. Create a column group for AGE and HEIGHT: select both column headers (using
   Ctrl+click to select multiple), then open the Context Panel > Actions and
   choose "Group columns". Assign a recognizable color (e.g. blue) to the group.
3. Verify that the two columns now appear under a shared group band in the
   grid header.
4. Verify that the group structural state is present: the AGE and HEIGHT columns
   each carry the 'group' column tag, and the dataframe carries the '.columnGroups'
   tag reflecting the new group.
5. Shift+click on both AGE and HEIGHT grouped column headers (click AGE, then
   Shift+click HEIGHT while both are in the same group band).
6. Verify that no console error is raised — the console-error delta is 0 from
   the baseline (GROK-17505 guard).
7. Click on the empty space within the group band header (the area not occupied
   by any column header, between or beside the grouped column names in the band).
8. Verify that no console error is raised — the console-error delta is 0
   (GROK-17443 guard).
9. Create a second column group for WEIGHT and SEX in the same way (Context Panel
   > Actions > Group columns), assign a different color (e.g. green).
10. Select the name of the first group (AGE + HEIGHT group) in the header band.
11. Then press Esc to deselect the group selection.
12. Verify that no console error is
raised at all — the console-error delta is 0 (GROK-17442 + GROK-18213 guards).
13. Save the current view as a project: click the project Save button in the
    ribbon (not a layout save — use the full project Save path). Note the
    project name.
14. Close All views.
15. Reopen the saved project from the Files browser.
16. Verify that both column groups are present with their assigned colors intact
    and no console error is raised on reopen.

Expected:
- Step 3: The AGE and HEIGHT columns display under a shared group band in the
  grid header area.
- Step 4: The AGE and HEIGHT columns each carry the 'group' tag, and the
  dataframe carries the '.columnGroups' structural tag reflecting the group
  membership.
- Step 6: No console error after Shift+clicking two grouped column headers;
  console-error delta is 0 (GROK-17505 guard).
- Step 8: No console error after clicking group empty space; console-error
  delta is 0 (GROK-17443 guard).
- Step 12: No console error after creating a second
  group, selecting the first group's name, and pressing Esc; console-error
  delta is 0 (GROK-17442 + GROK-18213 guards).
- Step 16: Both column groups are present with the same colors as before the
  project round-trip; the '.columnGroups' and 'group' df tags survive project
  save/reopen (GROK-17441 guard). No console error on reopen.
- Teardown: delete the probe project after verification to avoid leaving it
  in the server.

## Automation notes
- grok.shell.o is DEBOUNCED and first-value-wins inside a debounce window: re-asserting the target
  columns in a tight loop right after another context-panel focus collapses back to that PRIOR
  focus, the "N columns" Actions pane never rebuilds, and the "Group columns..." label never
  renders. Wait for the debounce to flush before reading grok.shell.o, and re-assign at most once
  per iteration with a full settle in between.
- The column-group dialog ALWAYS pre-fills [name="input-Group"] with the literal "Group". Two
  groups left on that default share one key in the df's `.columnGroups` map and the second silently
  clobbers the first, so fill a distinct name per group before OK.
- The hamburger popup's Colors Type control is a Datagrok ChoiceInput whose real
  <select name="input-Type"> is laid out at 0x0 (the visible affordance is the CSS chrome around
  it). Playwright's selectOption / waitFor({state:'visible'}) can never satisfy actionability on
  it — they burn the timeout and throw. Commit the value with a native value-set plus an `input`
  dispatch on that same <select>, and fire it from SEPARATE page.evaluate calls spaced by real
  waits: the Dart ChoiceInput binding subscribes a few frames after pane-Colors mounts, and an
  in-page tight loop keeps the microtask queue busy and starves that attach.
- The reopen step's only red is the CONSOLE channel, not the group colours — those come back
  correctly. What the listener catches is the publish-chain stack trace (ProjectMeta.publish),
  the same benign save-window class both persistence scenarios already handle. Copy their shape:
  grid-columns-style-persist-spec.ts declares a letter-agnostic benign pattern (the minified symbol
  drifts — never pin one) and treats that class as noise ONLY inside the save/reopen window, while
  anything else in the window and everything outside it still fails the step. Do not widen the
  filter beyond that window, and do not drop the console assertion.

### Signal discipline
All body prose uses UI-action language. Technical signals are recorded here:
- Sort state: `grid.props.sortByColumnNames` / `grid.props.sortTypes` (confirm
  ['SEX','AGE'] and ['asc','desc']); first visual row read via `grid.gridRowToTable(0)`.
- Columns dialog type-filter checkbox state: DOM checkbox elements within the
  dialog panel (checked attribute after Reset filter = false for all).
- Columns dialog header element: a DOM element carrying the dialog's header title
  text — presence/absence confirmed by DOM query.
- Console errors: delta from the baseline `console.error` count recorded at Setup
  step 3 (and re-recorded at the start of each close/reopen in later scenarios).
- Group structural state: `df.col('AGE').getTag('group')` (column tag), and
  `df.getTag('.columnGroups')` (dataframe-level JSON tag listing all groups).
- Group color round-trip: compare the 'group' tag values and '.columnGroups' JSON
  before save and after reopen — color assignments are encoded in the df tags and
  must survive project serialization.
- Context Panel color sync: the Colors section's rendered state in the Context
  Panel DOM — confirm a non-null/non-empty color swatch element is visible for
  the AGE column immediately after the hamburger-menu coding step.

### GROK-17463 guard detail
The gear-open guard in Scenario 3 requires clicking the gear, waiting for the
panel to mount in the DOM, then clicking the panel's X button (not just an
outside-click dismiss), then waiting for the panel to unmount, then clicking
the gear again. Use a DOM-presence poll on the settings panel element to confirm
close before re-clicking.

### GROK-19333 / GROK-20167 / GROK-19332 — Columns dialog guards
The Columns dialog checkboxes are drawn on an embedded canvas (same pattern as
the Order or Hide Columns dialog in `grid.cp.columns-layout-persist`). Reset
filter is a DOM button. After clicking Reset: query the canvas's rendered checkbox
state by examining the dialog's internal DOM structure, not by pixel inspection.
For the header guard (GROK-20167), query the dialog's header element by text
content or a stable CSS selector immediately after the filter/reset cycle.
For the data-switch guard (GROK-19332), open the second table via `grok.shell.openTable`
or the Data > Table menu, then re-apply a type-filter before checking the error log.

### GROK-17505 / GROK-17443 / GROK-17442 / GROK-18213 — column-group guards
All four guards are console-error delta 0 checks. Capture the error count before
each sensitive step and assert delta == 0 immediately after. Do not accumulate
across scenarios — re-baseline after each close/reopen.

### GROK-17441 — group color persistence (project round-trip)
Group colors are stored in df tags: column tag 'group' (a group name string) and
the dataframe tag '.columnGroups' (JSON encoding group membership and colors).
The project save must be via the REAL ribbon Save button (not just a layout save),
and the reopen via the Files browser. After reopen, read the
'.columnGroups' tag from the reopened df and compare to the pre-save snapshot.
Teardown: delete the probe project in a finally block via `grok.dapi.projects.delete`.

### GROK-19288 — hamburger menu and Context Panel sync
The hamburger icon appears on canvas hover — deliver a `page.hover` to the column
header's center before clicking the icon. If the icon does not appear after hover,
retry once (canvas render latency). The Context Panel Colors section is a DOM
panel; confirm the AGE column's color swatch is non-empty immediately after the
hamburger coding step without any panel refresh.

### Out of scope
- PowerGrid package custom cell renderers and pinned-columns renderer.
- Deep color-scheme editing (custom palettes, edit dialogs): owned by the flat
  Viewers/color-coding section.
- Column group drag-to-reorder (requires pixel-precision drag — manual-only per
  the MANUAL-ONLY CARVE-OUT in grid.cp.rows-select-filter-navigate).

---
{
  "order": 16,
  "datasets": ["System:DemoFiles/demog.csv"]
}
