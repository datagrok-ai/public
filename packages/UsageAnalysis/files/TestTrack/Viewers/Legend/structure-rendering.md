---
feature: legend
realizes_atlas: []
realizes: []
realized_as:
  - structure-rendering-spec.ts
target_layer: playwright
coverage_type: regression
priority: p2
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/structure-rendering.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - scenario-2-step-5-verification-depth
scope_reductions: []
related_bugs:
  - GROK-19083
---

# Legend — molecule structure rendering

Verifies that when a legend is bound to a column of chemical structures (`Molecule` type), it renders actual molecule thumbnails — not raw SMILES text or category indices — across seven viewer types, and that Scatter plot markers can also render as molecule glyphs. Also checks that structure rendering survives layout and project reloads.

## Setup

1. Open SPGI (`System:DemoFiles/chem/SPGI.csv`) and wait for `Molecule` semantic-type detection on the `Core` and `Structure` columns
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

- Relies on SPGI's `Core` and `Structure` columns being auto-detected as `Molecule` type.
- Standard legend UI flows (visibility, position, color picker, save-layout dialog behavior) are covered once in `visibility-and-positioning.md` and not exercised directly here.
- Setting the legend column is only a setup step here, not a test of column-switch behavior itself — that's covered in `visibility-and-positioning.md`.
- Bug GROK-19083 (deselecting markers desyncs the legend) has a dedicated repro (`legend-grok-19083-spec.ts`); this scenario verifies the positive baseline — with markers selected, the legend renders molecule glyphs correctly.

---
{
  "order": 2,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
