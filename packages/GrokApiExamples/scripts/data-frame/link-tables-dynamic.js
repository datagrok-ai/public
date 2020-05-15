// Master-detail table linking feature with dynamic data loading (see also link-tables.js).

let tickers = DataFrame.fromCsv(
`ticker
AAPL
TSLA
MSFT`);
grok.addTableView(tickers);

let detailsView = null;

tickers.onCurrentRowChanged.subscribe((_) => {
    grok.loadDataFrame(`https://www.quandl.com/api/v1/datasets/WIKI/${tickers.currentRow.ticker}.csv`)
        .then(function(t) {
            if (detailsView === null)
                detailsView = grok.addTableView(t);
            else
                detailsView.dataFrame = t;
        });
});
