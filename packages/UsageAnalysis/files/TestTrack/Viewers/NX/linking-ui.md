---
feature: viewers
realizes_atlas: []   # untyped: nx cross-viewer integration chain
priority: p2
target_layer: manual-only
coverage_type: smoke
---

Chain note: scenarios 1–5 of this folder build on each other and create five
server-side projects (NxProject, NxProjectCalcColumns, NxProjectViewers,
NxProjectFormulaLegend, NxProjectFiltering). Cleanup for all five lives at the
end of filtering-ui.md (order 5); if a run aborts mid-chain, delete whichever
of these projects already exist before re-running.

1. Open SPGI, SPGI-linked1, SPGI-linked2
1. Go to Data > Link Tables
1. Set Tables to SPGI and SPGI-linked1
1. Set Link type to 'selection to selection'
1. Set Key columns to Id <> Concept Id
1. Click Link
1. Click New Link
1. Set Tables to SPGI-linked1 and SPGI-linked2
1. Set Link type to  'selection to filter' 
1. Set Key columns to: 
   1. Sample Name
   1. link column 1
   1. link column 2
   1. link column 3
1. Add a line chart and set properties as follows:
   1. Data:
      1. Table to SPGI-linked2
      1. Filter to `${link column 3}=="v ii" && ${link column 1} <30`
      1. Split to link column 2
      1. Overview to link column 1
   1. X: X to Value1
   1. Y: Y Axis Type to logarithmic
1. Go to SPGI
1. Select/deselect some/all rows. Verify the line chart shows only the rows passing the link filter — nothing disappears or lingers unexpectedly as the selection changes
1. Save the project with datasync as NxProject
1. Close All
---
{
  "order": 1,
  "datasets": ["System:DemoFiles/chem/SPGI.csv","System:AppData/ApiTests/datasets/SPGI-linked1.csv","System:AppData/ApiTests/datasets/SPGI-linked2.csv"]
}