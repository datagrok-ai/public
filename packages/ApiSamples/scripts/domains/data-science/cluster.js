// An example of using clustering.
//
// https://datagrok.ai/help/explore/cluster-data

let table = await grok.data.loadTable('https://public.datagrok.ai/demo/xclara.csv');
let clusteredTable = await grok.ml.cluster(table, ["V1", "V2"], 3)
let view = grok.shell.addTableView(clusteredTable);
view.scatterPlot({
  x: 'V1',
  y: 'V2',
  color: 'clusters',
});