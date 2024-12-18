// https://datagrok.ai/help/visualize/viewers/shape-map

let t = await grok.data.loadTable('https://public.datagrok.ai/demo//earnings-by-state.csv');

let view = grok.shell.addTableView(t);
view.shapeMap();
