// Master-detail table linking feature with dynamic data loading (see also link-tables.js).

let tickers = DG.DataFrame.fromCsv(
  `chickweight 
fish`);
grok.shell.addTableView(tickers);

let detailsView = null;

tickers.onCurrentRowChanged.subscribe((_) => {
  grok.data.loadTable(`https://calmcode.io/static/data/${tickers.currentRow.ticker}.csv`)
    .then(function (t) {
      if (detailsView === null)
        detailsView = grok.shell.addTableView(t);
      else
        detailsView.dataFrame = t;
    });
});
