---
feature: legend
target_layer: playwright
coverage_type: smoke
priority: p0
pyramid_layer: ui-smoke
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/visibility-and-positioning.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - scenario-4-step-1-sp-legend-toggle-r-one-js-api-exposure
  - scenario-10-step-5-a-few-viewers-subset-for-mini-legend
  - scenario-11-step-3-verification-depth
scope_reductions: []
related_bugs:
  - GROK-17438
  - GROK-17222
  - github-3132
  - GROK-17278
  - GROK-19083
  - GROK-19041
---

# Legend — visibility, position, and core interactions

Smoke-tests the standard legend UI that every viewer shares: showing/hiding, positioning (including Auto and the four corner positions), resizing via the splitter, filtering by clicking a legend item (including exclude-via-cross-icon), the per-category color picker, the "(no value)" swatch for empty values, mini-legend mode, and persistence of all of this across layout and project reloads. Other Legend scenarios in this section delegate these standard flows here instead of repeating them.

## Setup

1. Open SPGI (`System:DemoFiles/SPGI.csv`) — wait for table view to be ready
2. Add seven viewers via the JS API: **Scatter plot**, **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**
3. Arrange viewers in the default tiled grid layout (`tv.addViewer(...)` placement is acceptable for automation; manual run can dock into a 2×4 grid)
4. Set the categorical legend on each viewer to `Stereo Category`, using the viewer's own legend-source property:
   * Scatter plot → **Color**
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**

## Scenarios

### 1. Legend appears on every viewer with matching swatch colors

1. Verify legend is present on every viewer — DOM `[name="legend"]` count equals 7
2. Verify swatch colors on each legend match the on-canvas points / bars / cells of the host viewer

### 2. Legend redraws on column change

1. On the Scatter plot, change the legend source (**Color** property) to `Series`
2. Verify the legend redraws with the new categories from `Series`
3. Change the legend source back to `Stereo Category`
4. Verify the legend redraws again with the original categories

### 3. Legend region resize via splitter

1. On any viewer, drag the legend splitter handle to enlarge the legend region (or set `viewer.props.legendWidth` / `legendHeight`)
2. Verify the legend region visibly grows
3. Drag the splitter again to shrink the legend region
4. Verify the legend region visibly shrinks

### 4. Filter rows via legend item click and cross-icon exclusion

1. Ctrl+click the `R_ONE` category in the Scatter plot legend (or call `sp.legend.toggle('R_ONE')` if the JS API exposes it)
2. Verify only `R_ONE` rows remain visible in the dataframe
3. Click the **X** icon on the `R_ONE` swatch to exclude it from the filter instead
4. Verify `R_ONE` is now excluded and the remaining categories are visible

### 5. Per-category color picker — Cancel / OK / propagation

1. Hover over a legend swatch on the Scatter plot to reveal the color picker icon, click it
2. Verify the color picker dialog opens
3. Change the color, click **Cancel**
4. Verify the change is discarded and the category retains its original color
5. Reopen the picker, change the color, click **OK**
6. Verify the new color is committed
7. Verify the new color propagates to every other viewer using `Stereo Category` for legend (the column is the single source of truth)

### 6. `(no value)` swatch on a column with empty values

1. Switch the legend source on the Scatter plot to `Primary Series Name`
2. Verify the `(no value)` / empty-value swatch is present in the legend
3. Hover over the `(no value)` swatch and click the color picker icon
4. Change its color via the dialog and click **OK**
5. Verify the color is committed and the swatch reflects the new color

### 7. Layout persistence — column, custom colors, visibility

1. Save the layout
2. Re-apply the saved layout (allow ≥3 s settle for the legend rebuild)
3. Verify the legend column, custom colors, and visibility state all persist after the round-trip

### 8. `Visibility = Always` and `Position = Auto`

1. Right-click the legend on every viewer and set **Visibility** = `Always`, **Position** = `Auto` (or set `viewer.props.legendVisibility = 'Always'`, `legendPosition = 'Auto'`)
2. Verify legends remain visible on all viewers
3. Resize a viewer (drag a split bar, or set `viewer.root.style.width = '300px'`)
4. Verify with `Position = Auto`, the legend repositions to whichever side has the most free space (Right vs Top vs Left vs Bottom per aspect ratio)
5. Save the layout, then re-apply it (allow ≥3 s settle)
6. Verify `Visibility = Always` and `Position = Auto` both persist

### 9. `Visibility = Auto` — auto-hide on small viewer

1. Right-click the legend on every viewer and uncheck auto-positioning, set **Visibility** = `Auto`
2. Reduce a viewer's size below ~250 px (`viewer.root.style.width = '200px'`)
3. Verify with `Visibility = Auto`, the legend hides on the small viewer (falls back to the mini-icon)
4. Restore the viewer size to ≥400 px
5. Verify the legend reappears

### 10. Corner positions and mini-legend mode

1. Set `legendPosition` on a viewer to `LeftTop`
2. Verify the legend renders in the left-top corner
3. Repeat for `LeftBottom`, `RightTop`, `RightBottom`
4. Verify each corner position renders correctly
5. Enable mini-legend mode (`viewer.props.miniLegend = true`) on a deterministic subset of viewers
6. Verify the mini variant is more compact than the full legend
7. Save the layout, re-apply it
8. Verify corner position and mini-legend mode both persist

### 11. Project persistence round-trip

1. Save the project
2. Close the project, then reopen it
3. Verify positioning state (corner positions + mini-legend mode + Visibility / Position settings) survives the persistence round-trip
4. Close All to clean up

## Notes

- This is the section's UI smoke test — all other Legend scenarios delegate the standard visibility/position/color-picker/splitter flows here instead of repeating them.
- Several bugs found in this area have their own dedicated repro specs and are not re-verified here: `legend-grok-17438-spec.ts`, `legend-grok-17222-spec.ts`, `legend-github-3132-spec.ts`, `legend-grok-17278-spec.ts`, `legend-grok-19083-spec.ts`.
- The legend resize handle (Scenario 3) and the color picker dialog (Scenario 5) need real UI interaction and can't be substituted with JS API calls.

---
{
  "order": 1,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
