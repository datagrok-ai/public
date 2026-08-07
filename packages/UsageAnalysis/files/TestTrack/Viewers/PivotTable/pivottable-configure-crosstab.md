---
feature: pivottable
realizes_atlas:
  - pivottable.cp.configure-crosstab-values
  - pivottable.int.agg-types-track-agg-columns
  - pivottable.int.empty-aggregates-clear-pivot
realizes:
  - viewers.pivot-table
priority: p0
target_layer: playwright
coverage_type: smoke
realized_as:
  - pivottable-configure-crosstab-spec.ts
related_bugs:
  - id: github-2535
    status: fixed
  - id: GROK-19114
    status: fixed
  - id: GROK-16899
    status: fixed
  - id: GROK-16305
    status: fixed
  - id: GROK-17841
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      The Group by row shows a chip for DIS_POP; the Aggregate row shows a chip
      for avg(AGE); the Pivot row shows a chip for SEVERITY. The aggregated
      workspace table opened by ADD contains DIS_POP as the key column,
      SEVERITY-derived columns as headers, and avg(AGE) values that match an
      independent df.groupBy(['DIS_POP', 'SEVERITY']).avg('AGE') computation.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      The column picker opens without any console-error delta (console.error
      count before == console.error count after): GROK-19114 guard.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      After clicking the Group by + and selecting SEX, a chip labeled SEX
      appears in the Group by row alongside DIS_POP. The aggregated table in the
      workspace updates: it now groups by both DIS_POP and SEX, and the row
      count equals the number of distinct (DIS_POP, SEX) pairs in demog.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      The avg(AGE) chip's context menu is open. A check mark is visible on the
      currently active aggregation item (avg / Mean) inside the Aggregation
      group. No console errors have been written since the right-click
      (GROK-17841 guard).
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      The check mark has moved from the previous item to Sum inside the
      Aggregation group. The chip caption updates to read sum(AGE). The
      aggregated column in the workspace table (or confirmed via JS API on the
      pivot's aggregated dataframe) reflects sum values. The menu is still open
      (radio items do not close the menu: GROK-16899 guard, first pick).
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      The check mark has moved again within Aggregation to a second item (e.g.
      Median). The chip caption reflects the new aggregation. The menu is still
      open (GROK-16899 guard, second pick in the same open menu).
  - anchor: "Scenario 2 Step 8"
    expectation: >-
      The chip caption reflects the newly picked column (e.g. avg(HEIGHT) or
      sum(HEIGHT)). The Aggregation group is rebuilt to show exactly the
      aggregation types that HEIGHT's numeric type supports — no unsupported
      items. The menu is still open. The previously active aggregation that was
      unsupported by HEIGHT has been replaced by a type-appropriate default (I2:
      pivottable.int.agg-types-track-agg-columns).
  - anchor: "Scenario 2 Step 10"
    expectation: >-
      Remove others closes the menu. Only one chip remains in the Aggregate row.
      No console errors since the Remove others click.
  - anchor: "Scenario 3 Step 3"
    expectation: >-
      The Aggregate row is empty (no chips). The Pivot row is hidden from the
      viewer (I4: pivottable.int.empty-aggregates-clear-pivot). Any pivot-column
      chips that were visible are silently cleared — the configuration is lost,
      not hidden.
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      The Aggregate row shows one chip (avg(AGE)). The Pivot row reappears in
      the viewer (the Pivot row is shown only while an aggregation exists).
  - anchor: "Scenario 4 Step 3"
    expectation: >-
      The property panel for the pivot viewer is open and its
      groupByColumnNames, pivotColumnNames, aggregateColumnNames and
      aggregateAggTypes fields are visible.
  - anchor: "Scenario 4 Step 5"
    expectation: >-
      Driving the aggregate column-list editor dialog from the property panel
      (search AGE, uncheck the first row; search WEIGHT, check the first row)
      fires the panel change handler: the editor count updates and
      aggregateColumnNames is rewritten to contain WEIGHT and to no longer
      contain AGE. No console error is written during the gesture (GROK-16305
      guard).
  - anchor: "Scenario 4 Step 7"
    expectation: >-
      Driving the group-by column-list editor dialog from the property panel
      (search DIS_POP, uncheck the first row; search SEX, check the first row)
      fires the panel change handler: the editor count updates and
      groupByColumnNames is rewritten to contain SEX and to no longer contain
      DIS_POP. No console error is written during the gesture (GROK-16305
      guard).
  - anchor: "Scenario 5 Step 4"
    expectation: >-
      After the Refresh click the visible tag chips reset to the type-driven
      auto-configuration: the Group by and Pivot rows are cleared and the
      Aggregate row shows avg(AGE) and avg(HEIGHT).
  - anchor: "Scenario 5 Step 9"
    expectation: >-
      After restoring the saved layout, the pivot viewer is a member of the
      current view's viewer set, and its groupByColumnNames, pivotColumnNames,
      aggregateColumnNames and aggregateAggTypes match the values that were set
      before saving (including the non-default aggregation): github-2535 guard.
  - anchor: "Scenario 5 Step 13"
    expectation: >-
      After closing all views and reopening the project, the pivot viewer is
      present in the view and its configuration (groupByColumnNames,
      pivotColumnNames, aggregateColumnNames, aggregateAggTypes) matches the
      saved state. Non-default aggregation survives the project round-trip
      (github-2535 guard).
---

# Pivot Table — Configure cross-tab values

## Setup

1. Log in to Datagrok.
2. Close all open views.
3. Open the demog dataset via the File Browser (Files > App Files > Demo Files > demog.csv)
   and wait for the table view to appear.
4. Add the Pivot Table viewer from the Toolbox (Viewers section).
5. Wait for the tag-editor header to appear with three rows labeled Group by, Aggregate, and Pivot.

## Scenarios

### Scenario 1: Default configuration, aggregation output, and column picker regression

Steps:
1. Observe the Group by, Aggregate, and Pivot rows in the viewer's tag-editor header.
2. Click the ADD button in the counts area at the top of the pivot viewer to publish the
   aggregated table to the workspace.
3. Wait for the new table to appear as an open tab in the workspace.
4. Read the new table's column names and at least one cell value from the first data row;
   independently compute df.groupBy(['DIS_POP', 'SEVERITY']).avg('AGE') via the JS API on
   the source dataframe and compare the results.
5. Record the current console error count.
6. Click the + icon next to the Group by row to open the column picker popup.
7. Close the column picker popup without selecting a column (press Escape or click outside it).
8. Read the current console error count and compare to the count recorded in Step 5.
9. Click the + icon next to the Group by row to open its column picker.
10. Select SEX from the column picker; wait for the SEX chip to appear in the Group by row.
11. Close the workspace table opened in Step 3 by clicking its tab close icon.

### Scenario 2: Aggregate tag context menu — multi-pick in one open menu

Steps:
1. Ensure the Aggregate row shows avg(AGE) (from the default auto-configuration; close and
   re-add the viewer if it was left in a different state from Scenario 1).
2. Record the current console error count.
3. Right-click the avg(AGE) chip in the Aggregate row to open its context menu.
4. In the Aggregation group of the context menu, note which item carries the check mark.
5. Click the Sum item (or any item other than the currently checked one) in the Aggregation
   group — do NOT close the menu.
6. Click a second item in the Aggregation group (e.g. Median) — do NOT close the menu.
7. Navigate to the Column group inside the same still-open context menu; click on a numeric
   column other than AGE (e.g. HEIGHT).
8. Read the Aggregation group — note the items now offered.
9. Click a second column in the Column group (e.g. WEIGHT) — do NOT close the menu.
10. Click Remove others in the context menu to close the menu and leave only the current chip.
11. Read the current console error count and compare to the count recorded in Step 2.

### Scenario 3: Pivot row hidden when last aggregate is removed

Steps:
1. Ensure the Aggregate row has at least one chip (avg(AGE)) and the Pivot row is visible
   with SEVERITY (from the default configuration; re-add the viewer if needed).
2. Click the × icon on every chip in the Aggregate row until the row is empty.
3. Observe the state of the Pivot row and the content of the viewer's tag-editor area.
4. Click the + icon next to the Aggregate row and select AGE from the column picker.
5. Observe the state of the Pivot row.

### Scenario 4: Property panel configuration path (GROK-16305 guard)

Steps:
1. Ensure the viewer is showing the default auto-configuration (DIS_POP in Group by,
   avg(AGE) in Aggregate, SEVERITY in Pivot); re-add if needed.
2. Record the current console error count.
3. Select the pivot viewer (click its header) so the viewer's property panel shows the
   column-list editors.
4. Locate the aggregateColumnNames field in the property panel and click it to open the
   column-list editor dialog.
5. In the dialog, type AGE into the search field — the list re-sorts so AGE is the first row —
   and clear the check box on that first row; then replace the search text with WEIGHT, and set
   the check box on the first row. Close the dialog to apply.
6. Wait for the Aggregate row chip in the viewer to update to the WEIGHT measure.
7. Locate the groupByColumnNames field in the property panel, click it to open its column-list
   editor dialog, search DIS_POP and clear the check box on the first row, then search SEX and
   set the check box on the first row. Close the dialog to apply.
8. Wait for the Group by row chip in the viewer to update to SEX.
9. Read the current console error count and compare to the count from Step 2.
10. Close the property panel.

### Scenario 5: Refresh resets configuration; layout and project persistence (github-2535)

Steps:
1. Close all views and reopen demog; add the Pivot Table viewer with the default
   auto-configuration (DIS_POP in Group by, avg(AGE) in Aggregate, SEVERITY in Pivot).
2. Read and record the current configuration: groupByColumnNames, pivotColumnNames,
   aggregateColumnNames, aggregateAggTypes.
3. Click the Refresh icon (the circular arrow) in the pivot viewer's command bar.
4. Verify the visible configuration resets to the type-driven defaults: the Group by and
   Pivot rows are cleared and the Aggregate row shows avg(AGE) and avg(HEIGHT).
5. Now set a non-default aggregation: right-click the avg(AGE) chip, open the context menu,
   and change the aggregation to Median (or Sum) — confirm the chip caption updates.
6. Save the current layout via the ribbon Save button in the toolbar, using the name
   "pivot-crosstab-test" (or equivalent), and wait for confirmation.
7. Change the Group by column to SEX by clicking × on DIS_POP and then + on Group by to
   add SEX — verifying the chip updates.
8. Re-apply the saved layout "pivot-crosstab-test" from the gallery and wait for the pivot
   to update.
9. Read the configuration (groupByColumnNames, pivotColumnNames, aggregateColumnNames,
   aggregateAggTypes) from the restored pivot.
10. Save the current workspace state as a project via the ribbon Save → Save as Project,
    naming it "pivot-crosstab-project".
11. Close all open views.
12. Reopen the project "pivot-crosstab-project" from the Projects panel.
13. Wait for the pivot viewer to appear and read its configuration.
14. Delete the layout "pivot-crosstab-test" and the project "pivot-crosstab-project" via
    their respective management panels (cleanup).

## Automation notes

- The aggregated-values comparison (Scenario 1 Step 4) is implemented as a JS-API check:
  `await page.evaluate(() => grok.shell.tableByName('demog').groupBy(['DIS_POP', 'SEVERITY']).avg('AGE').aggregate())` compared against the workspace table's values.
- The console-error delta checks (Scenarios 1, 2, 4) use a `page.on('console')` collector that
  records console messages of type error (plus page errors), ignoring the harmless cloned-iframe
  noise; the count is read before and after the guarded action, and a positive delta is a FAIL.
- The chip selector must use caption-based selection (not `.d4-tag` index) because tag-row rebuilds invalidate any pre-recorded index; confirmed live — the chip caption span survives rebuilds.
- The configuration read-back (Scenarios 4 and 5) is via `page.evaluate(() => { const v = grok.shell.viewers.find(v => v.type === 'PivotTableViewer'); return { groupBy: v.props.groupByColumnNames, pivot: v.props.pivotColumnNames, agg: v.props.aggregateColumnNames, aggTypes: v.props.aggregateAggTypes }; })`.
- Layout save and restore uses the real ribbon Save button and gallery UI, not a JS-API shortcut — the persistence channel is the UI path, not the API path (the distinction is the github-2535 regression guard).
- Project save/close/reopen uses `loginToDatagrok` after the reopen to ensure the session is active.
- Setup uses the `loginToDatagrok`, `closeAllViews` and `addViewerByIcon` helpers.
