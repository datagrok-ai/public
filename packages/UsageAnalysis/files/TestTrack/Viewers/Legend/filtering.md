---
feature: legend
realizes_atlas: []
realizes: []
realized_as:
  - filtering-spec.ts
target_layer: playwright
coverage_type: regression
priority: p1
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/filtering.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - scenario-8-step-5-6-filteredselected-behavior-with-empty-selection
  - scenario-4-step-2-cross-viewer-behavior-of-in-viewer-filter
  - scenario-9-step-3-deselect-a-few-ambiguity
scope_reductions: []
related_bugs:
  - GROK-17222
---

# Legend — stays in sync with filtering

Verifies that a viewer's legend always reflects the currently visible rows, no matter how those rows got filtered: the Filter Panel (numerical, categorical, structure filters), a viewer's own in-viewer filter property, click-to-filter interactions, and the different Row Source modes (All / Filtered / FilteredSelected / Selected). Checked across seven viewer types.

## Setup

1. Open SPGI (`System:DemoFiles/SPGI.csv`) — wait for table view + semantic-type detection on `Core` / `Structure` columns
2. Add seven viewers via the JS API: **Scatter plot**, **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**
3. Set the categorical legend on each viewer to `Stereo Category`, using the viewer's own legend-source property:
   * Scatter plot → **Color**
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**
4. Open the **Filter Panel**

## Scenarios

### 1. Filter Panel — numerical, categorical, and structure filters

1. Apply a numerical filter `Average Mass > 400` via the Filter Panel
2. Verify ≈ 1588 / 3624 rows remain (within tolerance)
3. Verify every viewer's legend lists only categories present in the filtered subset
4. Apply a categorical filter on `Stereo Category` keeping only `R_ONE` and `S_UNKN`
5. Verify legend on every viewer now lists exactly two categories
6. Apply a structure filter on `Core` (sketch a substructure or pick the first row's structure)
7. Verify viewer + legend update to reflect the structure-filtered subset

### 2. Layout persistence — filter state survives reload

1. Save the layout
2. Re-apply the layout — wait at least 3 s after `loadLayout` for the filter panel to rebuild
3. Verify the filter state and legends from Scenario 1 survive the round-trip

### 3. Reset filters

1. Reset all filters via `df.filter.setAll(true)`
2. Verify legends on every viewer return to the full category set

### 4. In-viewer filter property

1. On the Scatter plot, set the **Filter** property to `${Stereo Category} in ["R_ONE", "S_UNKN"]`
2. Verify only the two categories remain in the Scatter plot legend (filter applies in-viewer; other viewers may differ depending on viewer-level filter behavior)

### 5. Compose Filter Panel + in-viewer filter

1. With the in-viewer filter on Scatter plot still active, apply an additional Filter Panel filter `Average Mass > 300`
2. Verify both filters compose
3. Verify the Scatter plot legend stays at the two categories from the in-viewer filter

### 6. Click-to-filter on Bar / Pie / Trellis (and zoom-filter on Scatter)

1. On the Scatter plot, alt-drag a rectangle to zoom-and-filter (or set `sp.props.filter` to an x/y range expression)
2. Verify the Scatter plot legend drops categories with no points in the zoomed region
3. On the Bar chart, set **OnClick** = `Filter`, then click any bar
4. Verify only that bar's category remains in the Bar chart legend
5. On the Pie chart, set **OnClick** = `Filter`, then click any slice
6. Verify only that slice's category remains in the Pie chart legend
7. On the Trellis plot, set **OnClick** = `Filter`, then click any cell
8. Verify the Trellis legend reflects only the clicked cell's category

### 7. Layout persistence — click-to-filter state survives

1. Save the layout
2. Re-apply the layout (allow ≥3 s settle)
3. Verify the filter state from Scenario 6 is preserved across the round-trip

### 8. Row Source modes — legend reflects each row source

1. On the Scatter plot, set **Row Source** to `All`
2. Note legend categories shown — should reflect all rows
3. Set **Row Source** to `Filtered`
4. Verify legend reflects only filtered rows
5. Set **Row Source** to `FilteredSelected`
6. Verify legend reflects only filtered + selected rows
7. Set **Row Source** to `Selected`
8. Verify legend reflects only selected rows

### 9. Bar chart edge case — Stack with `includeNulls = false` [coverage_type: edge]

(Continues from Scenario 8 state on the Bar chart.)

1. On the Bar chart set: **Value** (count) = `CAST Idea ID`, **Category** = `Stereo Category`, **Stack** = `Primary Scaffold Name`
2. Uncheck **Value > Include nulls** (`bc.props.includeNulls = false`)
3. In the Filter Panel, deselect a few `Primary Scaffold Name` categories
4. Verify the Bar chart legend lists ONLY the stack categories still drawn — no ghost entries from the deselected ones

## Notes

- Covers the full filter-source matrix (Filter Panel x in-viewer filter x click-to-filter) across 7 viewer types. The Bar chart edge case (Scenario 9) additionally checks that unchecking "Include nulls" doesn't leave ghost legend entries for deselected stack categories.
- Standard legend UI flows (visibility, position, color picker, save-dialog behavior) are covered once in `visibility-and-positioning.md` and not repeated here.
- Scenario 9 continues from the Bar chart state left by Scenario 8 rather than starting fresh.
- Bug GROK-17222 ("Line chart: legend is not consistent with filtering") is the primary related bug — its dedicated repro (`legend-grok-17222-spec.ts`) reproduces the case where the legend misses a filter event. This scenario verifies the positive baseline: the legend should update on every filter source.

---
{
  "order": 3,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
