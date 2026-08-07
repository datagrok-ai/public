---
feature: pivottable
realizes_atlas:
  - pivottable.cp.chrome-history-and-drag-config
  - pivottable.int.history-menu-requires-existing-columns
  - pivottable.int.default-aggr-type-remembered
realizes:
  - viewers.pivot-table
priority: p2
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-17122
    status: fixed
  - id: GROK-16201
    status: fixed
  - id: GROK-16074
    status: fixed
  - id: github-3414
    status: fixed
  - id: GROK-14995
    status: fixed
realized_as:
  - pivot-table-spec.ts
ui_companion: pivot-table-ui.md
scope_reductions:
  - id: SR-01
    check: E-TRACE-02
    rationale: |
      Scenario 7's drag gestures are manual-only: the grid column headers are
      canvas (no DOM handle), the pivot uses d4 initDragDrop pointer tracking
      with no HTML5 drag hooks and no public drag/drop API, and the remove drop
      zone exists only during a live native drag — covered in pivot-table-ui.md.
    verdict_status: SCOPE_REDUCTION
expected_results:
  - anchor: "Scenario 1 Step 2"
    expectation: >-
      On a fresh viewer the auto-configuration reads back as Group by DIS_POP
      (the widest categorical column), Pivot SEVERITY (the second), Aggregate
      avg(AGE) (the first numerical column) — both as tag captions and as the
      groupByColumnNames / pivotColumnNames / aggregateColumnNames /
      aggregateAggTypes property values.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      After closing the viewer with the cross icon, the pivot is gone from the
      view's viewer set and the console error count is unchanged from the
      pre-close baseline (GROK-17122).
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      With Show Header set to false, the Data row container, the Group by /
      Aggregate / Pivot tag rows and the counts element are absent from the DOM
      (hidden for real, not merely restyled), and they all come back when the
      property returns to true.
  - anchor: "Scenario 3 Step 6"
    expectation: >-
      With Show Command Bar set to false, the command bar with the history and
      refresh icons disappears from the viewer, and reappears when the property
      returns to true.
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      The title set through the property panel is displayed in the viewer
      header, the description text appears above the grid when its position is
      Top, and it disappears when Description Visibility Mode is Never — each
      state read from the live DOM, with no substitution fallback.
  - anchor: "Scenario 5 Step 4"
    expectation: >-
      After Save parameters, localStorage['grok-aggregation-history'] holds a
      JSON array with an entry describing the saved configuration.
  - anchor: "Scenario 5 Step 6"
    expectation: >-
      Re-applying the saved history entry restores Group by, Aggregate and Pivot
      to the saved configuration, read back from the tag captions.
  - anchor: "Scenario 5 Step 8"
    expectation: >-
      After the WEIGHT column is removed from the table, the history menu no
      longer offers the entry whose configuration names WEIGHT (I8 —
      history-menu-requires-existing-columns).
  - anchor: "Scenario 6 Step 3"
    expectation: >-
      The Aggregate + picker reopens with no error; the remembered type renders
      on the canvas picker and is checked manually (I9 —
      default-aggr-type-remembered).
  - anchor: "Scenario 8 Step 2"
    expectation: >-
      Grouping by the identifier-like USUBJID column produces one aggregate row
      per identifier with no console error (GROK-16201).
  - anchor: "Scenario 8 Step 5"
    expectation: >-
      The aggregated result opened in the workspace keeps the semantic types of
      its key columns (GROK-16074).
  - anchor: "Scenario 8 Step 7"
    expectation: >-
      Switching the Table property value back and forth duplicates neither the
      table entries in the Data row (github-3414) nor the tag-editor header
      (GROK-14995).
---

# Pivot table — chrome, history and drag-driven configuration

## Purpose

Covers the pivot viewer's frame: what the viewer looks like when it appears,
how its chrome shows and hides, and the regressions pinned by GROK-17122, GROK-16201,
GROK-16074, github-3414 and GROK-14995. 
## Setup

1. Log in to Datagrok.
2. Close all open views.
3. Open the demog dataset via the File Browser (Files > App Files > Demo Files >
   demog.csv) and wait for the table view to appear.
4. Add the Pivot Table viewer from the Toolbox (Viewers section).
5. Wait for the tag-editor header with the Group by, Aggregate and Pivot rows.

## Scenarios

### Scenario 1: Default auto-configuration

Steps:
1. Look at the tag rows of the freshly added viewer.
2. Verify the auto-configuration: Group by shows DIS_POP, Pivot shows SEVERITY,
   Aggregate shows avg(AGE); the same three lists read back from the viewer's
   properties (groupByColumnNames, pivotColumnNames, aggregateColumnNames with
   aggregateAggTypes).
3. Verify the counts area is visible below the tag rows.

### Scenario 2: Add and remove the viewer (GROK-17122)

Steps:
1. Record the current console error count.
2. Close the pivot viewer with the cross icon on its title bar.
3. Verify the viewer is gone from the view's viewer set and the console error
   count equals the baseline from Step 1.
4. Re-add the Pivot Table viewer from the Toolbox and verify the same
   auto-configuration as in Scenario 1 comes back.

### Scenario 3: Show Header and Show Command Bar

Steps:
1. Open the viewer's property panel.
2. Set Show Header to false.
3. Verify the Data row, the Group by / Aggregate / Pivot tag rows and the
   counts area disappear from the viewer.
4. Set Show Header back to true and verify they reappear.
5. Set Show Command Bar to false.
6. Verify the command bar with the history and refresh icons disappears.
7. Set Show Command Bar back to true and verify it reappears.

### Scenario 4: Title and description display modes

Steps:
1. In the property panel, set Show Title to true and Title to "My Pivot".
2. Verify "My Pivot" appears in the viewer header.
3. Set Description to "Summary stats" and Description Position to Top.
4. Verify the description text appears above the pivot grid.
5. Set Description Visibility Mode to Never and verify the description
   disappears.
6. Restore Show Title and the description properties to their defaults.

### Scenario 5: Aggregation history

Steps:
1. Configure Group by RACE, Aggregate avg(WEIGHT), no Pivot.
2. Click the history icon in the command bar and choose Save parameters.
3. Reconfigure to Group by SEX, Aggregate avg(AGE).
4. Verify localStorage['grok-aggregation-history'] contains a JSON array with
   an entry for the RACE / avg(WEIGHT) configuration.
5. Open the history menu and pick the saved RACE / avg(WEIGHT) entry.
6. Verify Group by returns to RACE and Aggregate to avg(WEIGHT).
7. In the main grid, remove the WEIGHT column (column context menu > Remove),
   then close the pivot viewer and re-add it from the Toolbox — the history
   menu re-filters its entries when the viewer is re-attached.
8. Open the history menu again and verify the RACE / avg(WEIGHT) entry is no
   longer offered.
9. Click the refresh icon in the command bar and verify Group by and Pivot are
   cleared and Aggregate shows avg(AGE) and avg(HEIGHT). Re-open demog to
   restore the WEIGHT column before the next scenario.

### Scenario 6: Remembered aggregation type

Steps:
1. Click + next to Aggregate, pick HEIGHT and choose sum in the aggregation
   popup — a sum(HEIGHT) tag appears.
2. Remove the sum(HEIGHT) tag with its cross icon.
3. Click + next to Aggregate again and verify the popup pre-offers sum as the
   aggregation type before any choice is made.
4. Cancel the popup and verify the Aggregate row is unchanged.


### Scenario 8: ID grouping, semantic types and the Table property

Steps:
1. Set Group by to USUBJID (the subject identifier column), Aggregate to
   avg(AGE), no Pivot.
2. Verify the inner grid shows one aggregate row per subject identifier and no
   console error is written.
3. Set Group by back to DIS_POP.
4. Click ADD in the counts area to open the aggregated result in the
   workspace.
5. Verify the key column of the opened table keeps its semantic type from the
   source column.
6. Close the aggregated table.
7. In the Data row, switch the Table property to another value and back, and
   verify the Data row shows each table entry once (no duplicates) and the
   tag-editor header is not duplicated.
8. Close all open views to clean up.

## Automation notes

- Console-error baseline and delta: count console messages of level error
  before and after the action; the cloned-iframe warning class is excluded.
- localStorage read: `await page.evaluate(() =>
  window.localStorage.getItem('grok-aggregation-history'))`, parsed as JSON.
- The viewer-set membership read:
  `Array.from(grok.shell.tv.viewers).map(v => v.type)`.
- Semantic type read-back: `grok.shell.t.col('<name>').semType` on the opened
  aggregated table.
- Confirmed live: the counts caption starts with the aggregated row count
  ("N rows"); history-menu entries are captioned by the aggregation set (e.g.
  "key(RACE),avg(WEIGHT)"). The inner-grid drag leg is manual — see
  pivot-table-ui.md.
- The Aggregate `+` popup and the tag context menus are canvas/flyout
  structures — drive them with the sustained-hover pattern documented in
  grok-browser/references/viewers/pivot_table.md and realized in
  pivottable-configure-crosstab-spec.ts.
---
{
  "order": 17,
  "datasets": ["System:DemoFiles/demog.csv"]
}
