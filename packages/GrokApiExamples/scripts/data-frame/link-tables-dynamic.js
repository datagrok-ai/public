// Master-detail table linking feature with dynamic data loading (see also link-tables.js).

let tickers = DataFrame.fromCsv(
`ticker
AAPL
TSLA
MSFT`);
gr.addTableView(tickers);

var detailsView;

tickers.onCurrentRowChanged(function (_) {
    gr.loadDataFrame(`https://www.quandl.com/api/v1/datasets/WIKI/${tickers.currentRow.ticker}.csv`)
        .then(function(t) {
           if (detailsView == null)
               detailsView = gr.addTableView(t);
           else
               detailsView.dataFrame = t;
        });
});