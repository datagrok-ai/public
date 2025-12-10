const df = await grok.data.files.openTable('System:AppData/ClinicalCase/ae.csv');
let view = grok.shell.addTableView(df);
view.addViewer('timelines');