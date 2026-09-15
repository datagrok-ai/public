const df = grok.data.demo.demog();
const tv = grok.shell.addTableView(df);
tv.grid.columns.byName('subj').pin();
