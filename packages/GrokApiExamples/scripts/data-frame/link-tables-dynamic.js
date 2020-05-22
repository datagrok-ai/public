// Master-detail table linking feature with dynamic data loading (see also link-tables.js).

let tickers = DG.DataFrame.fromCsv(
`ticker
AAPL
TSLA
MSFT`);
grok.shell.addTableView(tickers);

let detailsView = null;

tickers.onCurrentRowChanged.subscribe((_) => {
    grok.data.loadDataFrame(`https://www.quandl.com/api/v1/datasets/WIKI/${tickers.currentRow.ticker}.csv`)
        .then(function(t) {
            if (detailsView === null)
                detailsView = grok.shell.addTableView(t);
            else
                detailsView.dataFrame = t;
        });
});
