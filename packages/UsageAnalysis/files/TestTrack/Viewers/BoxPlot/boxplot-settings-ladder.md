---
feature: boxplot
realizes_atlas:
  - boxplot.cp.value-category-axes-persist
  - boxplot.int.category-sets-marker-color
realizes:
  - viewers.box-plot
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-18876
    status: fixed
  - id: GROK-18515
    status: fixed
  - id: GROK-19752
    status: fixed
  - id: GROK-18456
    status: fixed
  - id: GROK-20469
    status: fixed
  - id: GROK-20395
    status: fixed
  - id: GROK-20397
    status: fixed
realized_as:
  - boxplot-settings-ladder-spec.ts
expected_results:
  - anchor: PRE-LADDER
    expectation: >-
      With Value set to STARTED (a datetime column), the Axis Type control in
      the Context Panel > Value section is disabled or visually inactive (the
      dependsOn gate for a non-numerical/datetime value is active; GROK-20395).
      The Axis Type menu item in the context menu stays present and enabled —
      the gate is property-panel-only.
  - anchor: Step 4
    expectation: >-
      After setting Category 1 to SEX, the Marker Color Column in Context Panel
      > Color automatically switches to SEX with categorical coloring — the
      category-sets-marker-color side effect fires
      (boxplot.int.category-sets-marker-color).
  - anchor: Step 5
    expectation: >-
      After explicitly setting Marker Color Column to HEIGHT and then changing
      Value to WEIGHT, the Marker Color Column remains HEIGHT and the
      invertColorScheme, colorMin, and colorMax settings are untouched — an
      axis-column change must not override an explicitly chosen coloring column
      (GROK-18876).
  - anchor: Step 8
    expectation: >-
      After setting Axis Type to Log, no error or warning appears in the browser
      console (GROK-18515 — the log-axis switch must not throw). The viewer's
      visible range stays within the positive data bounds with no stretched
      empty space beyond the data range (GROK-20397).
  - anchor: Step 10
    expectation: >-
      After zooming with the range slider so the viewport narrows and then
      setting Marker Color Column to a different column, the range slider
      position and the viewer viewport are unchanged — setting marker color must
      not reset the range slider (GROK-20469).
  - anchor: Step 11
    expectation: >-
      After the group-comparison tail (showGroupComparison on, control group
      chosen, Adjust By covariate set), the Context Panel shows the Adjust by
      caption and a non-empty covariate selector.
  - anchor: LAYOUT-ROUND-TRIP
    expectation: >-
      After re-applying the saved layout, the Box Plot viewer is present in the
      layout and the Scatter Plot that was added as the re-arm step is absent.
      Every setting from the ladder (Value=WEIGHT, Category 1=SEX, Category
      2=RACE, Show Minor Categories on, Show All Categories on, Marker Color
      Column=SEX (set in Scenario 2 Step 18, after the explicit HEIGHT phase),
      invertColorScheme on, colorMin=20, colorMax=80, Axis Type=Log, Invert Y
      Axis on, Plot Style=violin) reads back from the restored viewer's Context
      Panel. The viewport zoom is NOT asserted on layout restore (viewport is
      excluded from layout by design).
  - anchor: PROJECT-ROUND-TRIP
    expectation: >-
      After reopening the saved project (Close All, then reopen), the Box Plot
      viewer is present with the full ladder configuration AND the viewport zoom
      is restored to match the state before save (GROK-19752, GROK-18456 — zoom
      state must survive the project round-trip). showGroupComparison,
      controlGroup, and covariateColumnName read back to the values set in the
      group-comparison tail.
---

# Box Plot — Settings Ladder and Persistence Round-Trip

## Setup

1. Close all open tables and viewers.
2. Open the demog golden dataset: navigate to Files > App Data and open
   System:DemoFiles/demog.csv.
3. Add a Box Plot viewer to the table view via the toolbar (Add viewer
   > Box Plot).

## Scenarios

### Scenario 1: Pre-ladder datetime gate probe (GROK-20395)

Steps:
1. In Context Panel > Value, set Value to STARTED (a datetime column).
2. In Context Panel > Value, look at the Axis Type control.
3. Verify that the Axis Type control is disabled or visually inactive in the
   property panel — Axis Type is unavailable for datetime value columns
   (GROK-20395). Do NOT verify the context-menu Axis
   Type item; that item stays present and enabled in the menu — the gate
   is property-panel-only.
4. In Context Panel > Value, set Value back to AGE to begin the ladder.

Expected:
- With a datetime Value column selected, the Axis Type control in the Context
  Panel is disabled (the property-panel gate for non-numerical/datetime values
  is active; GROK-20395). The context-menu Axis Type item is unaffected.

### Scenario 2: Settings ladder — accumulate, do not revert

Steps:
1. In Context Panel > Value, confirm Value is AGE (set it if needed).
2. In Context Panel > Data, set Category 1 to SEX.
3. Verify that the Marker Color Column in Context Panel > Color automatically
   changed to SEX with categorical coloring (setting a category auto-selects
   the marker color while no explicit coloring is active).
4. In Context Panel > Color, set Marker Color Column explicitly to HEIGHT
   (overrides the auto-selected SEX).
5. In Context Panel > Color, enable Invert Color Scheme, set Color Min to 20,
   and set Color Max to 80.
6. In Context Panel > Data, set Category 2 to RACE (two-level X axis).
7. In Context Panel > Data, enable Show Minor Categories.
8. In Context Panel > Data, enable Show All Categories.
9. In Context Panel > Value, change Value to WEIGHT.
10. Verify that the Marker Color Column is STILL HEIGHT and that Invert Color
    Scheme, Color Min, and Color Max are untouched — an axis-column change must
    not steal the explicitly chosen coloring column (GROK-18876).
11. In Context Panel > Value, set Value Min to 20 and Value Max to 60.
12. In Context Panel > Value, set Axis Type to Log.
13. Verify that no error or warning appears in the browser console after the
    Axis Type switch to Log (GROK-18515 — the log-axis switch must not throw)
    AND that the viewer's visible range stays within the positive data bounds
    with no stretched empty space above the data (GROK-20397).
14. In Context Panel > Value, enable Invert Y Axis.
15. In Context Panel > Style, set Plot Style to violin.
16. Use the range slider at the bottom of the Box Plot to zoom in (drag the
    slider handles inward to narrow the visible value range).
17. Verify that the viewer viewport changed (the visible value range is narrower
    than before the drag).
18. In Context Panel > Color, change Marker Color Column to SEX (a different
    column from HEIGHT, while the zoom is in place).
19. Verify that the range slider position and the viewer viewport are unchanged
    after setting the Marker Color Column — marker color must not reset the range
    slider (GROK-20469).

Expected:
- After setting Category 1 to SEX, the Marker Color Column automatically
  switches to SEX (category-sets-marker-color interaction fires).
- After explicitly setting Marker Color Column to HEIGHT and then changing Value
  to WEIGHT, the Marker Color Column remains HEIGHT and the color-scheme settings
  (Invert Color Scheme, Color Min, Color Max) are untouched (GROK-18876).
- After setting Axis Type to Log, no console error or warning appears and the
  viewer's visible range stays within the positive data bounds (GROK-18515,
  GROK-20397).
- After zooming with the range slider and then changing Marker Color Column, the
  range slider position and viewport are unchanged (GROK-20469).

### Scenario 3: Group-comparison persistence tail

Steps:
1. In Context Panel > Data, enable Show Group Comparison.
2. In the on-chart control strip, select a control group from the control group
   selector (pick any available category in the control selector).
3. In Context Panel > Group Comparison, set Adjust By to HEIGHT (covariate
   column).
4. Verify that the "Adjust by:" caption and a non-empty covariate selector are
   visible in the Context Panel, confirming the group-comparison state is set.

Expected:
- After enabling Show Group Comparison, choosing a control group via the
  on-chart control selector, and setting an Adjust By covariate, the Context
  Panel shows the "Adjust by:" caption and a non-empty covariate selector —
  the group-comparison state is set and will be persisted.

### Scenario 4: Layout round-trip (without viewport)

Steps:
1. Save the current view layout using the ribbon View > Save Layout button or
   via the Layout menu (note the layout name for later retrieval).
2. To re-arm the layout (modify the current layout away from the saved state):
   close the Box Plot viewer and add a Scatter Plot viewer to the same table
   view.
3. Confirm the current view now shows a Scatter Plot (no Box Plot).
4. Re-apply the saved layout using the ribbon View > Apply Layout menu or the
   Layouts panel (select the layout saved in Step 1).
5. Verify that the Box Plot viewer is present in the layout AND the Scatter Plot
   that was added as the re-arm step is absent (the layout restores the exact
   viewer set from the save point).
6. Verify that every setting from the ladder reads back from the restored Box
   Plot's Context Panel: Value=WEIGHT, Category 1=SEX, Category 2=RACE, Show
   Minor Categories on, Show All Categories on, Marker Color Column=SEX (the
   value set in the group-comparison tail's preceding step), Invert Color Scheme
   on, Color Min=20, Color Max=80, Axis Type=Log, Invert Y Axis on, Plot
   Style=violin.
7. Delete the probe layout from the Layouts panel so it does not persist across
   runs (even on failure, clean up the layout in a finally-equivalent step).

Expected:
- After re-applying the saved layout, the Box Plot viewer is present and the
  Scatter Plot added as the re-arm step is absent.
- Every setting from the ladder reads back from the restored viewer's Context
  Panel: Value, Category columns, Marker Color Column, color-scheme settings,
  Axis Type, Invert Y Axis, Plot Style.
- The viewport zoom is NOT required to be present on layout restore — the
  viewport is excluded from layouts by design.

### Scenario 5: Project round-trip (with viewport and group-comparison props)

Steps:
1. Save the current project using the ribbon's Save button (Save > name
   the project > confirm).
2. Close All (close the demog table and all viewers).
3. Reopen the saved project from the Projects panel or via File > Projects.
4. Verify that the Box Plot viewer is present with the full ladder configuration
   — read every setting from the Context Panel as in the layout round-trip, AND
   additionally verify the viewport zoom is restored (the visible value range
   matches the zoomed state from before the save; GROK-19752, GROK-18456).
5. Verify that the group-comparison persistence props are restored: Show Group
   Comparison is on, the control group is still selected, and Adjust By shows
   HEIGHT (the covariate is restored).
6. Delete the probe project from the Projects panel (even on failure, clean up
   in a finally-equivalent step).

Expected:
- After Close All and reopening the saved project, the Box Plot viewer is
  present with the full ladder configuration — all settings from the ladder
  read back as saved.
- The viewport zoom is restored to match the pre-save zoomed state (GROK-19752,
  GROK-18456 — zoom state must survive the project round-trip).
- Show Group Comparison is on, the control group is still selected, and Adjust
  By shows HEIGHT after project reopen (group-comparison persistence).

## Automation notes

- target_layer rationale: the settings ladder requires driving Context Panel
  property inputs across multiple categories (Value, Data, Color, Style,
  Group Comparison); the layout and project round-trips require UI actuation
  of the Save / Apply Layout and Save Project ribbon buttons; Playwright is
  the only layer that can perform all of these in a live Datagrok session.
- Viewport scope: bp.viewport is @Prop(includeInLayout: false)
  (box_plot_look.dart:247-248). Assert viewport ONLY in the project round-trip
  (Scenario 5). Do NOT assert the viewport on layout re-apply (Scenario 4).
- category-sets-marker-color side effect fires only while no explicit user
  coloring is active; after an explicit Marker Color Column set it is
  suppressed for subsequent category changes.
- Log-axis signals: grok.shell.warnings is not exposed to JS — the browser
  console/pageerror delta across the Axis Type=Log switch is the no-throw
  channel (GROK-18515); the GROK-20397 bound check reads the value-axis
  viewport against the positive data max.
- Project save must go through the real ribbon Save button
  (helpers/projects.ts saveProjectViaUI), not the JS API.
