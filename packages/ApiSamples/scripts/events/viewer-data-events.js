//name: viewer-data-events
//language: javascript
//description: Handling viewer events: onDataEvent, onDataSelected, onDataHovered, onDataRowClicked.
let view = grok.shell.addTableView(grok.data.demo.demog());
let viewers = [view.scatterPlot(), view.histogram(), view.pieChart(), view.barChart(), view.pcPlot()];

for (let v of viewers) {
  v.onDataEvent.subscribe((e) =>
    grok.shell.info(e.type + ': ' + (e.row ? (e.row + ' row') : (e.bitset.trueCount + ' rows'))));
}

grok.shell.info('Hover or click data to see event details. Different viewers provide different details.');