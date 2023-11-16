const df = await grok.data.files.openTable('Samples:Files/energy_uk.csv');
let view = grok.shell.addTableView(df);
view.addViewer('chord');