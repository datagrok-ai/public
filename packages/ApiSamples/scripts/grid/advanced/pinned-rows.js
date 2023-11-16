const df = grok.data.demo.demog();
const tv = grok.shell.addTableView(df);

// setting the column-value pairs
tv.grid.setOptions({
  pinnedRowValues: ['1', '0', '4'],
  pinnedRowColumnNames: ['subj', 'subj', 'subj'],
});