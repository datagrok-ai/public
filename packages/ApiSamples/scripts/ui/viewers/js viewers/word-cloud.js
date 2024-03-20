const df = await grok.data.files.openTable('System:DemoFiles/word_cloud.csv');

let view = grok.shell.addTableView(df);

view.addViewer('word cloud');