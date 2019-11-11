gr.loadDataFrame('https://www.quandl.com/api/v1/datasets/WIKI/AAPL.csv')
    .then(t => gr.addTableView(t));
