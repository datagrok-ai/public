// Once specified, custom category order gets used for sorting across all visualizations

let t1 = DG.DataFrame.fromCsv(
    `Severity,id
Medium, s_001
High, s_002
High, s_003
Medium, s_004
Low, s_005`);

t1.columns.byName('severity').setCategoryOrder(['Low', 'Medium', 'High']);

let view = grok.shell.addTableView(t1);

// Grid will have its values sorted appropriately
view.grid.sort(['severity']);

// Values on the "Severity" axis are appropriately sorted
view.scatterPlot({title: 'Custom sort order'});
