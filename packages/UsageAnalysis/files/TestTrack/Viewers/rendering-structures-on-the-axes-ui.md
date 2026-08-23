---
feature: viewers
realizes_atlas: []   # untyped: cross-viewer axis rendering, no single-feature atlas
priority: p2
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Rendered molecule images on canvas axes require visual judgment.
---

# Rendering structures on the axes — manual checklist

## Structures on axes and range sliders ([#2330](https://github.com/datagrok-ai/public/issues/2330))

1. Close all, open SPGI
2. For each of the following viewers — Scatter plot, Box plot, Bar chart,
   Trellis plot, Pie chart:
   - set the category axis (X / category / split, as applicable) to Structure
   - verify molecule structures render as axis/category labels, not empty
     labels or raw SMILES text
3. On each of the five viewers, drag the axis range sliders — structures
   re-render for the visible range, the axis never goes blank, no errors in
   the console (regression for #2330)

## Cell renderer handlers (Med Chem)

1. Close all, open Med Chem.csv (or the Med Chem project)
2. Click cells with custom renderers — for each, the Context Panel shows the
   info matching that cell's semantic type
3. In the Jira column, enter `GROK-11111` — the cell renders the ticket
   summary and additional info
4. Click that cell — the Context Panel shows the ticket details and the link
   to the ticket opens correctly

---
{
  "order": 15,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
