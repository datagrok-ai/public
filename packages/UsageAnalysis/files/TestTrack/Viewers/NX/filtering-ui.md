---
feature: viewers
realizes_atlas: []   # untyped: nx cross-viewer integration chain
priority: p2
target_layer: manual-only
coverage_type: smoke
---

1. Open the NxProjectFormulaLegend project (saved by the previous scenario) - it opens without errors
1. Go to the first SPGI view
1. Open the Filter Panel and add a scaffold tree filter by Structure
1. Upload the old scaffold tree `scaffold-tree-for-testing.tree` (QA-owned fixture: a scaffold tree saved from a pre-NX build; ask the QA owner for the file) - the scaffolds are not colored
1. Upload another old scaffold tree `scaffold-tree-for-nx-testing.tree` (QA-owned fixture: a scaffold tree saved from an NX build) - the scaffolds are colored
1. Check some checkboxes, change another setting on the scaffold tree - check the filtering (including other tabs) 
1. Save the layout
1. Close the Filter Panel - the filtering by the scaffold tree should be reset
1. Apply the saved layout - the filtering should be applied
1. Go to another SPGI view
1. Right-click any molecule in the Structure Column and select Current Value > Use as Filter - check the Filter Panel and the filtering (including the other tabs)
1. View > Layout > Clone View
1. On the new view, change the structure in the filter - check the filter on the previous view
1. Click the `core` column’s header
1. Go to the Context Panel > Chemistry > Rendering
1. Set Filter Type to `categorical`
1. Open the Filter Panel - the Core filter should be categorical 
1. View > Layout > Clone View
1. Set up filtering so that some rows are filtered out
1. Turn off all filters - all filtering should be turned off, all rows displayed.
1. Turn on the filters again
1. Turn off some individual filter - rows filtered out by this individual filter should be shown again
1. Apply the layout saved at step 7, and additionally a layout saved by an earlier NX scenario - the filtering each of them carries is applied, no errors
1. Save a new copy of the project as 'NxProjectFiltering'
1. Close All
1. Open the NxProjectFiltering project - it opens with the saved filtering applied, no errors
1. Cleanup (end of the NX chain): delete the five projects created by scenarios 1-5 - NxProject, NxProjectCalcColumns, NxProjectViewers, NxProjectFormulaLegend, NxProjectFiltering. After a partial run, delete whichever of them exist.
---
{
  "order": 5,
  "datasets": ["System:DemoFiles/chem/SPGI.csv", "System:AppData/ApiTests/datasets/SPGI-linked1.csv", "System:AppData/ApiTests/datasets/SPGI-linked2.csv"]
}
