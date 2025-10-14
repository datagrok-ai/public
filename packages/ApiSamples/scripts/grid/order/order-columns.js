// Order columns. Columns not mentioned in the parameters will be positioned after the specified ones.

let view = grok.shell.addTableView(grok.data.demo.demog());
view.grid.columns.setOrder(['age', 'sex', 'race']);
