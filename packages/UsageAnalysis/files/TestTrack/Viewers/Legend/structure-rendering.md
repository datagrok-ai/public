---
feature: legend
sub_features_covered:
  - legend.show-main-item-icons
  - legend.column
  - legend.item.marker-picker
  - legend.refresh.on-data-change
target_layer: playwright
coverage_type: regression
priority: p1
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/structure-rendering.md
migration_date: 2026-05-07
migration_report: structure-rendering-migration-report.md
related_bugs:
  - GROK-19083
---

## Setup

1. Open SPGI (`System:DemoFiles/SPGI.csv`) and wait for `Molecule` semantic-type detection on the `Core` and `Structure` columns
2. Add seven viewers via the JS API: **Scatter plot**, **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**

## Scenarios

### 1. Molecule thumbnails render in legend across 7 viewer types

1. On each viewer, set the legend column to `Core` (a `Molecule` semType column), using the viewer's own legend-source property:
   * Scatter plot → **Color**
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**
2. Verify every legend renders rendered molecule thumbnails (not raw SMILES text or category indices)

### 2. Marker glyphs render as molecules on Scatter plot

1. On the Scatter plot, set **Marker** = `Core` and **Color** = `Core`
2. Verify markers on the Scatter plot canvas render as molecule glyphs (not default shape markers)
3. Change **Color** to `Series` (keep **Marker** = `Core`)
4. Verify markers stay rendered as molecule glyphs
5. Verify the legend now shows `Series` colors alongside the structure markers

### 3. Layout persistence — molecule thumbnails survive reload

1. Save the layout
2. Re-apply the saved layout (allow ≥3 s settle)
3. Verify every legend still renders structures after the round-trip (no regression to SMILES text or indices)

### 4. Project persistence round-trip

1. Save the project
2. Close the project, then reopen it
3. Verify legends and structure rendering survive the persistence round-trip (molecule thumbnails on every legend, molecule glyphs on Scatter plot markers)

## Notes

- Specialty: molecule thumbnail + marker glyph rendering across 7 viewers. Cross-cutting on `Molecule` semantic-type detection in `Core` / `Structure` columns of SPGI.
- Delegates standard legend UI flows (visibility / position / color picker / save-layout dialog widgets) to `visibility-and-positioning.md`. This scenario does NOT exercise those flows directly.
- `pyramid_layer: integration` — multi-subsystem (7 viewer types × molecule semantic type rendering). Not a smoke; not source-matrix; not bug-focused.
- `legend.show-main-item-icons` is the primary atlas sub-feature; `legend.column` is touched as a setup step (not actively tested for column-switch behavior — that lives in the smoke).
- Bug GROK-19083 (markers deselect ↔ legend sync) is cross-cutting; the bug-focused spec `legend-grok-19083-spec.ts` (proposed by chain) reproduces the deselect-markers invariant; this scenario verifies the positive baseline (markers selected → legend renders glyphs).

---
{
  "order": 2,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
