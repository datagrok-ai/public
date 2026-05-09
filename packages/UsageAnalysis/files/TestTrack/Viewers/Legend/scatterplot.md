---
feature: legend
sub_features_covered:
  - legend.column
  - legend.extra-column
  - legend.show-main-item-icons
  - legend.item.color-picker
  - legend.allow-item-coloring
  - legend.refresh.on-data-change
  - legend.item.click
  - legend.color-scale.numerical
  - legend.use-custom-color-coding
target_layer: playwright
coverage_type: regression
priority: p1
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/scatterplot.md
migration_date: 2026-05-07
migration_report: scatterplot-migration-report.md
related_bugs:
  - GROK-17438
  - GROK-17222
  - github-3132
  - GROK-17278
  - GROK-19083
---

## Setup

Each Scenario is independent (starts with "Open SPGI" and ends with "Close All" cleanup) — see Notes for spec strategy.

## Scenarios

### 1. Color + Marker combined legend on Scatter plot [coverage_type: edge]

1. Open SPGI (`System:DemoFiles/SPGI.csv`)
2. Add a scatterplot
3. Set **Color** = `Series` and **Marker** = `Series`
4. Verify the legend is combined (one legend entry per `Series` category showing both color and marker glyph)
5. Verify the color picker icon appears on hover; change colors for a few categories
6. Save the layout, then re-apply it
7. Verify color changes persist after the layout round-trip
8. Save the project, close it, and reopen it
9. Verify color changes persist across the project round-trip
10. On the plot, add a new **Color** column with a categorical formula:
    `if(${Stereo Category}=='S_UNKN', null, ${Series})` — verify a categorical legend appears
11. Set **Color** = `ID` and **Marker** = `Core` — verify the legend is displayed correctly
12. Close All

(Visual check «linear (numerical) legend with gradient swatch on `if(${Stereo Category}=='S_UNKN', null, ${Average Mass})`» moved to `legend-ui.md` §4, 2026-05-08 — when a numerical Color column is bound, the Scatter plot's `[name="legend"]` element is empty in the DOM and the gradient swatch renders to canvas without DOM children.)

### 2. Legend updates on X-axis change with derived nullable columns

1. Open SPGI
2. Add new columns:
   * `col1`: `if(${Stereo Category}!='S_UNKN', null, ${Average Mass})`
   * `col2`: `if(${Stereo Category}=='S_UNKN', null, ${Average Mass})`
3. Add a scatterplot
4. Set **X axis** = `col1`, **Color** = `Stereo Category`
5. Verify the legend reflects categories present in `col1`-non-null rows
6. Change **X axis** = `col2`
7. Verify legend categories update according to the new data subset (different `Stereo Category` distribution)
8. Test zooming/filtering on the new X-axis
9. Verify legend stays consistent with the visible data subset
10. Close All

### 3. In-viewer filtering — multiple scatterplots with shared filter [coverage_type: edge]

1. Open SPGI
2. Add a scatterplot
3. Set **Marker** = `Stereo Category`
4. Apply in-viewer filter `${Stereo Category} in ["R_ONE", "S_UNKN"]`
5. Verify the legend reflects only the two filtered categories
6. Add a second scatterplot with the same in-viewer filter and same Marker setting
7. Save the layout, then re-apply it
8. Verify both scatterplots' legends still reflect the filtered subset after the round-trip
9. Close All

### 4. Filter Panel filtering and click-to-filter on Scatter plot legend

1. Open SPGI
2. Add a scatterplot
3. Set **X/Y axes** = `Chemical space X` / `Chemical space Y`
4. Set **Color** = `Primary scaffold name`, **Marker** = `Stereo Category`
5. In the Filter Panel, apply a filter on `Primary scaffold name` (deselect a few categories)
6. Verify the data and legend update — legend lists only the categories still drawn
7. Click `R_ONE` in the scatterplot legend
8. Verify additional filtering is applied — no previously-filtered-out points reappear (legend filter composes with Filter Panel filter, not replaces)
9. Close All

### 5. Color coding from grid — linear and categorical, with persistence [coverage_type: edge]

1. Open SPGI
2. Add a scatterplot, a box plot, and a PC plot
3. Set **Color** = `Chemical Space X` on the scatterplot
4. In the grid, enable **linear color coding** for `Chemical Space X`
5. Verify the scatterplot and box plot legends reflect the linear color scheme (numerical color scale, not categorical items)
6. Change the color schema (pick a different palette), then invert it, then apply it to the text
7. Verify legends update after each schema change
8. Save the layout, then re-apply it
9. Verify the color changes persist after the layout round-trip
10. Change the grid's `Chemical Space X` coding from linear to **categorical**, then modify a few colors
11. Verify legends update — scatterplot now shows categorical legend items; box plot follows suit
12. Save the project, close it, and reopen it
13. Verify the categorical color customizations persist across the project round-trip

## Notes

- 5 independent sub-sections. Strategy `chained_tests` per chain — single spec file with 5 `test()` blocks (one per Scenario), each with its own `Open SPGI` setup and `Close All` cleanup. No state shared between Scenarios.
- Specialty: scatterplot-specific legend behaviors (Color/Marker combined, derived nullable color columns, in-viewer filter, click-to-filter, grid color coding linear-vs-categorical, persistence).
- Delegates standard legend UI flows (visibility / position / save-dialog widgets) to `visibility-and-positioning.md`.
- `pyramid_layer: integration` per chain Rule 4 — multi-subsystem (scatterplot × multiple legend behavior dimensions).
- 5 cross-cutting bugs intersect this scenario:
  - GROK-17438 (color picker → legend hides on shared-legend viewers) — Scenario 1 step 5.
  - GROK-17222 (legend not consistent with filtering) — Scenario 4 step 5.
  - github-3132 (sequential color changes reset previous) — Scenario 1 step 5.
  - GROK-17278 (linechart colors not saved to layout/project) — Scenario 1 steps 6–9 (positive baseline; linechart-specific in line-chart.md).
  - GROK-19083 (markers deselect ↔ legend sync) — Scenario 1 step 11 (Color=ID, Marker=Core); Scenario 3 step 3 (Marker=Stereo Category).
- Helpers: `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login`. No new helpers proposed.
- Visual / pixel-level check split into companion `legend-ui.md` (`target_layer: ui-only`, manual QA) on 2026-05-08: §4 «numerical color gradient swatch in scatter legend» (Sc 1 step 10 first sub-bullet). Spec body retains JS-API state assertions (`colorCodingType: 'Linear'`, `colorColumnName` round-trip) for the related flows.

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
