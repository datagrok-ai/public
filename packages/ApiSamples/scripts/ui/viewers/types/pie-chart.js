// https://datagrok.ai/help/visualize/viewers/pie-chart

let view = grok.shell.addTableView(grok.data.demo.demog());

view.pieChart();
view.pieChart({category: 'race', mode: 'Donut', centerLabel: 'Race'});
