
const df = await grok.data.files.openTable("Demo:Files/energy_uk.csv");
let view = grok.shell.addTableView(df);
view.addViewer('sankey');