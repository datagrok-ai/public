---
feature: viewers
realizes_atlas: []   # untyped: nx cross-viewer integration chain
priority: p2
target_layer: manual-only
coverage_type: smoke
---

1. Open the NxProject project
1. Go to SPGI
1. View > Layout > Clone View
1. Add new column:
   * `if(${Whole blood assay 1} != null, ${Whole blood assay 1}, if(${Route Admin}=="PO", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}=="Dog", 30.9, if(${Species}=="Monkey", 43.6, if(${Species}=="Minipig", 39, null)))))*100` - **there should be a warning about if(num,qnum)**
   * Enter another formula: `if(${NIBR logP} != null, ${NIBR logP}, if(${Route Admin}=="PO", ${Whole blood assay 1} / ${Chemical Space X} * 100 / 6 / ${Average Mass} * 1000000.0,null))/if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, if(Contains(${Species}, 'Mouse'), 125, if(${Species}=="Dog", 30.9, if(${Species}=="Monkey", 43.6, if(${Species}=="Minipig", 39, null)))))*100`
   * Name: `${Species} result`
1. Order Or Hide Columns: 
   * Unselect all
   * Select columns: NIBR logP, Whole blood assay 1, Chemical Space X, Average Mass, Species, and  `${Species} result`
1. Go to the `${Species}` result hamburger menu and click Add Filter - **the Filter Panel should open with the `${Species}` result filter**
1. Go to the filter tab and click missing values > filter out missing values
1. Rename the `Species` column to `Spec` - the `${Species} result` column should be renamed to `${Spec} result`
1. Change some values in visible columns - the corresponding value of the `${Spec} result` column should be changed
1. View > Layout > Clone View - the view with exactly the same set of visible columns and filter should open
1. Right-click the grid and select Add > Summary columns > Pie chart
1. Save the Layout
1. Change the layout (e.g. close the Filter Panel and resize the grid)
1. Apply the saved layout - check the filtering, `${Spec} result`, and summary column
1. Right-click the view’s header and select Table > Add view - a new view without opened Filter Panel and with all SPGI columns visible should appear; check: new grid should be filtered to 30 rows as all views are.
1. Open the Filter Panel and use Remove All
1. Use the Select Columns.. dialog to add Structure and some categorical filters - the Structure filter is the last in the filter list
1. Drag the Structure filter tab to the top of the Filter Panel 
1. Reset the filter  - filtering should be reset on all views
1. Save a new copy of the project as 'NxProjectCalcColumns'

Close All

Cleanup: deferred to the end of the NX chain (see filtering-ui.md).
---
{
  "order": 2,
  "datasets": ["System:DemoFiles/chem/SPGI.csv", "System:AppData/ApiTests/datasets/SPGI-linked1.csv", "System:AppData/ApiTests/datasets/SPGI-linked2.csv"]
}
