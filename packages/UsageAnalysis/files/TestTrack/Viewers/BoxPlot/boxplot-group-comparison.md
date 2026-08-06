---
feature: boxplot
realizes_atlas:
  - boxplot.cp.group-comparison-ladder
  - boxplot.cp.covariate-adjust-baseline
  - boxplot.int.covariate-sets-adjustment
realizes:
  - viewers.box-plot
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - boxplot-group-comparison-spec.ts
related_bugs: []
expected_results:
  - anchor: "Scenario 1 Step 1"
    expectation: >-
      The pre-comparison baseline holds: showPValue reads true and
      showGroupComparison reads false (int prop reads — the plain t-test
      overlay is the starting state).
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      After clicking the reveal icon (name='show-group-stats'), the viewer prop
      showGroupComparison reads true — a DOM-actuated prop state transition
      confirmed via int prop read.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      Hovering the overall p-value element yields a DOM tooltip whose text
      carries the test conclusion (e.g. 'Welch's t-test.').
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After switching the method via the on-chart select inside
      .d4-box-plot-group-comparison-controls, the prop echo confirms the new
      method value (secondary signal; primary is the DOM actuation).
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      After switching Category 1 to RACE (3+ groups), the group-comparison
      overlay recomputes to one-way ANOVA — the settle-gated canvas repaint IS
      the signal (waitForCanvasChange over a captured baseline; the prop echo
      is logging only).
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      After picking a control group in the on-chart control select, the prop
      controlComparisons reads true — int signal. The control band and per-group
      p row are present as canvas floors (no-error / settle-gated diff).
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      A new DataFrame appears in grok.shell.tables whose name contains the test
      prefix (e.g. 't-Test:') AND contains 'by RACE' as two separate substring
      checks. The table carries expected column structure and p-values are read
      from it. The context panel shows the BoxPlotStatsMeta accordion on
      per-group p row click (DOM).
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      The group-comparison overlay shows THREE effect rows for the Category1 ×
      Category2 interaction (canvas floor). A new DataFrame appears in
      grok.shell.tables whose name contains 'Two-Way ANOVA:' AND contains 'by
      RACE, SEX' as two separate substring checks; effect rows are read from the
      table.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      After clicking name='close-group-stats', all group-comparison UI is hidden
      AND the plain t-test p overlay returns (transition pair —
      showGroupComparison is false, plain showPValue overlay visible).
  - anchor: "Scenario 2 Step 2"
    expectation: >-
      After setting covariateColumnName=HEIGHT, the prop adjustmentMode reads
      'regressOut' (int side-effect signal). The DOM shows the 'Adjust by:'
      caption with HEIGHT and a 'Regress-out' toggle. The value-axis title
      remains the raw column name (no 'WEIGHT (adj)' or 'WEIGHT / HEIGHT'
      suffix).
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      After flipping the on-axis toggle to Ratio, adjustmentMode reads 'ratio'
      (int) and the DOM toggle reads 'Ratio'.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      After setting method=ANCOVA via the on-chart select, the prop method reads
      'ANCOVA' (uppercase literal, int). The Ratio/Regress-out toggle is hidden
      (zero-size / absent in the DOM) while the 'Adjust by:' caption stays
      visible. The stored adjustmentMode is unchanged.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      A new DataFrame appears in grok.shell.tables whose name contains the
      ANCOVA prefix AND 'by RACE' as two separate substring checks (the adjusted
      value may insert an infix so value+by adjacency is never asserted;
      equality on the full name is never asserted); ANCOVA results are read from
      the table columns.
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      With Category2 set and baseline switched to 'Matched · per stratum',
      per-stratum control bands appear as a canvas floor. The Simpson's-paradox
      red '!' cue appears deterministically next to the baseline toggle on the
      engineered fixture; its hover tooltip carries the Simpson message (DOM).
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      After clearing the covariate (covariateColumnName=''), the props
      adjustmentMode reads '' AND the method is cleared (no ANCOVA) — confirmed
      via int reads (box_plot_core.dart:534-537 master).
---

# Box Plot — Group Comparison and Covariate Adjustment

Acceptance path for the group-comparison ladder (GROK-20173) and covariate
adjustment (GROK-20381), both Ready for QA. These are new features, not
bugs — they are deliberately absent from related_bugs.

## Setup

1. Close all open views.
2. Open the demog dataset (System:DemoFiles/demog.csv).
3. Add a Box Plot viewer to the table view.
4. In the Context Panel > Value, set Value to AGE.
5. In the Context Panel > Data, set Category 1 to SEX (2 groups).
6. Ensure the viewer is sized sufficiently for the statistics table to be
   visible (the table auto-hides entirely when values do not fit — size
   the viewer before any statistics or comparison step).

## Scenarios

### Scenario 1: Group-comparison ladder

The group-comparison reveal icon appears only when the row count is below 200k,
the number of category combinations is between 2 and 30, and a value column
with data is present. With SEX (2 groups) and AGE on demog, these conditions
are satisfied.

Steps:

1. Verify the plain t-test p-value overlay is visible on the viewer canvas
   (the p-value overlay is on by default with 2 or more categories; this is the
   pre-comparison baseline — no tooltip-text check here, see Step 4 for the
   tooltip).
2. Hover the bare p-value area on the canvas.
3. When the reveal icon appears at the top of the viewer, click it to enable
   group comparison. Verify the group-comparison strip becomes visible on the
   viewer (confirming the feature activated).
4. Hover the overall p element in the group-comparison strip. Verify a tooltip
   appears carrying the test conclusion (e.g. "Welch's t-test." for 2 groups).
5. In the group-comparison controls bar on the chart, locate the method
   selector and switch the method from Welch to Student. The recomputed p value
   is verified from the results table in Step 8.
6. In the Context Panel > Data, set Category 1 to RACE (3+ groups, switches
   to one-way ANOVA). Confirm the group-comparison strip updates to the ANOVA
   method (visible as a change on the chart).
7. On the control group selector visible on the chart (visible on hover until a
   control is chosen), pick a control group (e.g. 'Caucasian'). Use the
   on-chart selector — setting the control group via the Context Panel alone
   does not activate control comparison. Verify the control band and per-group
   p row appear on the chart.
8. Right-click the control band or stats region and select "Add Control
   Comparisons Table" from the context menu. Verify:
   - A new table appears in the workspace.
   - Its name contains the test prefix (e.g. 't-Test:' or 'Control
     Comparisons:') and also contains 'by RACE' as separate parts of the name
     (do not assert these two parts appear directly adjacent — an adjusted value
     may appear between them; do not assert the full name exactly).
   - The table contains columns for group, p, and similar statistics.
   - Read p-values from the table as the primary numerical check.
   Click the per-group p row in the stats strip. Verify the Context Panel
   updates to show the statistics details accordion.
9. In the Context Panel > Data, set Category 2 to SEX (two-way ANOVA).
   Wait for the recompute — the overlay may briefly show the stale one-way
   result after the category is set; wait for the viewer to settle before
   checking (toggling group comparison off and back on is the recovery path).
   Verify the overlay shows THREE effect rows for RACE, SEX, and RACE × SEX.
   Right-click and select "Add Two-Way ANOVA Table". Verify:
   - A new table appears in the workspace.
   - Its name contains 'Two-Way ANOVA:' and also contains 'by RACE, SEX' as
     separate parts.
   - Effect rows are read from the table (primary check).
10. Click the close icon at the top of the viewer to dismiss group comparison.
    Verify:
    - All group-comparison UI is hidden (the controls bar and stats overlay are
      gone from the viewer).
    - The plain t-test p overlay returns (the transition back to the default
      p-value display).

### Scenario 2: Covariate adjustment and control-baseline modes

Start fresh from the Setup state (Value=AGE, Category1=SEX). This scenario
builds a separate ladder — do NOT carry any group-comparison state from
Scenario 1.

Steps:

1. In the Context Panel > Value, change Value to WEIGHT. In the Context
   Panel > Data, set Category 1 to SEX. Enable group comparison by hovering
   the bare p-value area and clicking the reveal icon that appears. Confirm the
   group-comparison strip is active.
   Important: set the covariate (Step 2) BEFORE picking ANCOVA as the method
   (Step 4). The regress-out mode activates only when the covariate is set
   while the method is not yet ANCOVA.
2. In the Context Panel > Data (or via the covariate selector on the chart),
   set Adjust By to HEIGHT. Verify:
   - The 'Adjust by:' caption with HEIGHT appears persistently on the chart.
   - The Ratio/Regress-out toggle is present and reads 'Regress-out'.
   - The value-axis title stays the raw column name (e.g. 'WEIGHT') —
     do NOT expect 'WEIGHT (adj)' or 'WEIGHT / HEIGHT'; the adjusted column
     is indicated only by the caption and mode toggle, not by the axis label.
3. Flip the Ratio/Regress-out toggle on the chart to Ratio. Verify:
   - The toggle reads 'Ratio'.
4. In the Context Panel > Data, set Category 1 to RACE (3+ groups). Pick a
   control group in the on-chart control group selector (e.g. 'Caucasian'). In
   the group-comparison controls bar on the chart, set the method to ANCOVA.
   Verify:
   - The method shown in the controls bar reads 'ANCOVA'.
   - The Ratio/Regress-out toggle is hidden while the 'Adjust by:' caption
     remains visible.
   - The adjustment mode set in Step 3 is unchanged (ANCOVA suppresses its
     effect but does NOT clear the stored value).
   - Do NOT verify the reverse direction (setting adjustment mode while ANCOVA
     is active resetting the method) — that path is not reachable through the
     UI when ANCOVA is selected.
5. Right-click the stats region and select "Add ANCOVA Table". Verify:
   - A new table appears in the workspace.
   - Its name contains the ANCOVA test prefix and also contains 'by RACE' as
     separate parts (do not assert these parts appear directly adjacent or assert
     the full name exactly).
   - The table contains ANCOVA result columns (F, df, p, adjusted values) and
     values are read from those columns (primary check).
6. In the Context Panel > Data, set Category 2 to SEX. Switch the baseline
   toggle on the chart to 'Matched · per stratum'. Verify per-stratum control
   bands appear on the chart. Simpson's-paradox fixture: create two calculated
   columns via Add New Calculated Column with engineered OPPOSITE within-stratum
   trends (remove them at the end of the test). Verify the red '!' cue appears
   next to the baseline toggle and its hover tooltip carries the
   Simpson's-paradox message. If the fixture cannot be driven reliably, note the
   result as inconclusive.
7. Clear the covariate: set Adjust By to None. Verify:
   - The 'Adjust by:' caption and mode toggle are no longer visible on the
     chart.
   - The method is no longer ANCOVA.

## Automation notes

- The 'Add ... Table' name contract is built by _tableNameSuffix
  (box_plot_group_comparison.dart:1823 master) and
  box_plot_stats_presentation.dart:271-415 master. Assert two separate
  `.contains()` checks — one for the test prefix and one for 'by
  <categories>' — never a full-string equality.
- statsInfoAt is implemented at box_plot_group_comparison.dart:1912 master;
  the context panel accordion is the BoxPlotStatsMeta component.
- Simpson's-paradox cue: box_plot_group_comparison.dart:43 (threshold) and
  :123 (cue emission) master.
- The 'Adjust by:' caption and mode toggle visibility is driven by
  _showAdjust (box_plot_group_comparison.dart:794 master) — hidden when
  ANCOVA is active.
- The adjustment-clears-ANCOVA code (box_plot_group_comparison.dart:681-683
  master) lives ONLY in the on-chart toggle's handler, which is hidden while
  ANCOVA is active — UNREACHABLE via props or menu. Do NOT author steps
  asserting that direction.
- target_layer rationale: the group-comparison and covariate surfaces are
  driven through on-chart native selects (the group-comparison controls bar,
  CSS class .d4-box-plot-group-comparison-controls), hover-reveal icons
  (name='show-group-stats' / name='close-group-stats'), right-click context
  menus, and Context Panel props — all requiring a live browser; no
  callable entry point exists for these flows.
- The group-comparison reveal icon carries name='show-group-stats'; the close
  icon carries name='close-group-stats'. The on-chart method/control/baseline
  selects are DOM elements; hover tooltips and the statistics detail accordion
  in the Context Panel are DOM. The workspace table list is grok.shell.tables.
- Prop echo signals: showGroupComparison (bool), controlComparisons (bool),
  adjustmentMode ('regressOut' | 'ratio' | ''), method ('Student' | 'Welch' |
  'ANCOVA' | ''), covariateColumnName (string). ANCOVA must be uppercase literal
  — a lowercase 'ancova' resets to ''. Covariate clear uses
  covariateColumnName='' (box_plot_core.dart:534-537 master).
- Numerical results channel: the primary channel for p / F / df / adjusted p
  values is the "Add ... Table" workspace table — read column structure and
  values from it. Canvas-drawn p values, control bands, and per-group p rows
  are secondary no-error floors because canvas digits are not readable by
  automation.
- Control group comparison activates ONLY when a control group is chosen in
  the on-chart selector. Setting the control group via the Context Panel alone
  does NOT activate it. Always drive through the on-chart selector.
- Two-way ANOVA: after setting Category 2 the overlay may briefly show the
  stale one-way result — wait for the viewer to recompute before asserting the
  three-row structure.
- Add-Table aiming recipe: the exclusive strip menu with the
  'Add ... Table' item opens only when the right-click lands INSIDE the
  comparison strip, which starts ~40 CSS px right of the canvas origin —
  clicks at smaller offsets (e.g. +20,+20) open the FULL viewer menu, which
  has no Add items. Recompute the canvas rect fresh before every attempt,
  right-click at canvas origin +42,+16 (fallback offsets rightward/below,
  e.g. +50,+16 / +60,+14 / +80,+14 / +42,+30 all land on the strip), then
  click the menu item by its name= attribute (div-Add-t-Test-Table /
  div-Add-Two-Way-ANOVA-Table / div-Add-Control-Comparisons-Table /
  div-Add-ANCOVA-Table) — a synthetic element .click() is sufficient.
- The added results table opens its OWN TableView and makes it current —
  after the new table appears, switch back to the home TableView (find it in
  grok.shell.tableViews by dataFrame name) before any further box-plot step,
  or every later viewer/DOM lookup dangles on the results view.
- In the added tables the p-value column is a STRING column with formatted
  values ('.0038', '<.0001') — parse the numeric out of the string instead of
  expecting a double from col.get().
- Simpson's-paradox fixture recipe (demog.csv, Welch t):
  create two calculated columns via addNewCalculated —
  SIMPSON_STRAT = `if(Mod(Round(${HEIGHT} * 1000), 2) == 0, "A", "B")`
  (HEIGHT sub-mm digit parity: balanced ~50/50 and sex-independent) and
  SIMPSON_VAL = `if(${SEX} == "M", 0, if(Mod(Round(${HEIGHT} * 1000), 2) == 0,
  2.5, -2.5)) + (Mod(Round(${WEIGHT} * 137), 600) / 30)` (for F a +2.5/-2.5
  effect by stratum, for M zero, plus wide 0..20 noise decorrelated from sex
  by the mod-wrap). Bind Value=SIMPSON_VAL, Category1=SEX, Category2=
  SIMPSON_STRAT, control='M', baseline 'matched'. Reference stats: stratum A
  diff +2.7 (t=+9.8), stratum B diff -2.7 (t=-11.9), pooled F-vs-M t=-0.18 —
  opposite significant within-stratum arms with a canceled pooled effect,
  which is the product's Simpson trigger (per-group signs contain both +1 and
  -1, at least one arm significant, pooled non-significant). The red '!' cue
  is the only red span (color rgb(235, 103, 103)) in the viewer; it sits just
  left of the baseline toggle and is display-gated. WEIGHT decimal parity is
  NOT usable as a stratum source (many integer weights skew it ~70/30), and
  low-multiplier noise transforms (e.g. HEIGHT*731 mod 2000) keep a ~1.0
  sex bias that makes the pooled test significant — both rejected. Tear both
  columns down in finally and assert the removal.
