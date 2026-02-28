const df = await grok.data.files.openTable('System:AppData/ApiSamples/ae.csv');
let view = grok.shell.addTableView(df);
view.addViewer('timelines');