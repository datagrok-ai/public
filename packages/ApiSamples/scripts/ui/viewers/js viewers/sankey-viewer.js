const df = await grok.data.files.openTable('System:DemoFiles/energy_uk.csv');
let view = grok.shell.addTableView(df);
view.addViewer('sankey');