---
feature: pcplot
realizes_atlas:
  - pcplot.cp.panel-chart-filter-composition
realizes:
  - viewers.pc-plot
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - pcplot-filter-panel-composition-spec.ts
related_bugs: []
expected_results:
  - anchor: "panel filter drops the count"
    expectation: "df.filter.trueCount drops below the full row count after a
      Filter Panel histogram filter (AGE 30-50) is applied"
  - anchor: "in-chart slider drops it further"
    expectation: "df.filter.trueCount drops below the Filter-Panel-only value
      after narrowing the HEIGHT axis range slider on the PC Plot — the two
      filter sources combine with AND"
  - anchor: "Reset View returns to the panel-only value"
    expectation: "df.filter.trueCount returns EXACTLY to the Filter-Panel-only
      value after Reset View on the PC Plot — Reset View clears ONLY the
      in-chart part and the panel filter survives"
  - anchor: "re-drag drops it again"
    expectation: "df.filter.trueCount drops below the Filter-Panel-only value
      once more after re-narrowing the HEIGHT range slider"
  - anchor: "Reset filters restores the full count"
    expectation: "df.filter.trueCount returns to the full row count after the
      Filter Panel Reset filters button — it clears BOTH filter sources"
---

# PC Plot — Filter Panel + In-Chart Filter AND-Composition and Reset Scoping

## Purpose

Verify that a Filter Panel filter and the PC Plot's own in-chart range slider
combine with AND on the shared DataFrame, and that the two reset paths have
different scope: Reset View on the PC Plot clears only the in-chart part (the
panel filter survives), while the Filter Panel Reset filters button clears
both sources. Every check reads the filtered row count the product computes,
never a slider setting read back (see Automation notes for how the count is
read).

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a PC Plot viewer to the current table view via the toolbar (Add viewer > PC Plot).
4. Set the axes to AGE, HEIGHT, WEIGHT (the axes carry range sliders).
5. Open the Filter Panel.
6. Record the full filtered row count (all rows pass the filter, read dynamically).

## Scenarios

### Scenario 1: AND-composition ladder with reset scoping

Steps:
1. Apply a Filter Panel histogram filter on AGE (window 30-50) — the filtered row count drops below the full count.
2. Narrow the HEIGHT range slider on the PC Plot (drag its max handle downward) — the filtered row count drops further: the Filter Panel and in-chart filters combine with AND.
3. Right-click the PC Plot > **Reset View** — only the in-chart filtering resets; the filtered row count returns exactly to the Filter-Panel-only value (the panel filter survives).
4. Narrow the HEIGHT range slider again — the filtered row count drops below the Filter-Panel-only value once more.
5. On the Filter Panel click **Reset filters** — both filter sources reset and the filtered row count returns to the full count.

Expected:
- The panel filter drops the count below the full row count
- The in-chart slider drops it further (AND-composition of the two sources)
- Reset View returns the count exactly to the panel-only value — the panel filter survives
- Re-dragging the slider drops the count again
- Reset filters restores the full count — both sources are cleared

## Automation notes

- "The filtered row count" is read dynamically as `df.filter.trueCount` on the
  shared dataframe, never hard-coded.
- The Filter Panel histogram filter is applied via
  `grok.shell.tv.getFiltersGroup().updateOrAdd({type: 'histogram', column: 'AGE', min: 30, max: 50})`;
  the Reset filters button is
  `.d4-filter-group-header [name="icon-arrow-rotate-left"]`.
- The in-chart filter is driven by dragging the real
  `[name="axis-slider-HEIGHT"] [name="max-handle"]` DOM element with standard
  mouse events.
- Reset View is actuated through the canvas context menu (`Reset View` item) —
  the scriptable-headless equivalent of the whitespace double-click; it fully
  restores the in-chart filter.

---
{
  "order": 8,
  "datasets": ["System:DemoFiles/demog.csv"]
}
