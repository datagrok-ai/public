// Hides completely empty columns

let tv = grok.shell.getTableView('myData');
for (let c of tv.dataFrame.columns)
  if (c.stats.missingValueCount === tv.dataFrame.rowCount)
    tv.grid.col(c.name).visible = false;