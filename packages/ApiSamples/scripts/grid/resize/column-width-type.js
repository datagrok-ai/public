// Apply a predefined column-width policy to the grid.
// Available options: Minimal, Compact, Optimal, Maximal.
// Equivalent to the right-click menu: Column Sizing.

const view = grok.shell.addTableView(grok.data.demo.demog());

view.grid.setColumnsWidthType(DG.ColumnWidthType.Maximal);

// Apply a different policy to a single column:
view.grid.columns.byName('age').setWidthType(DG.ColumnWidthType.Minimal);
