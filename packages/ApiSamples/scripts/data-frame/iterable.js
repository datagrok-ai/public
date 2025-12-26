// Uses ColumnList as iterable and prints names of columns

let demog = grok.data.demo.demog();

for (let column of demog.columns)
  grok.shell.info(column.name);