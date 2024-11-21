
const demog = grok.data.demo.demog(100);
demog.columns.setOrder(['height', 'weight', 'age']);
grok.shell.addTableView(demog);