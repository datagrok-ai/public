// The following demo shows, how to work with dynamically changing data in Datagrok.

function runStockBroker() {
  const count = 50;
  let tickers = {
    'AAPL': 200,
    'TSLA': 200,
    'GOOG': 200,
    'NFLX': 200
  };

  let df = DG.DataFrame.create(0);
  df.name = 'stocks';
  df.columns.addNew('symbol', DG.TYPE.STRING);
  df.columns.addNew('time', DG.TYPE.DATE_TIME);
  df.columns.addNew('price', DG.TYPE.FLOAT);
  df.columns.addNew('volume', DG.TYPE.FLOAT);
  df.columns.addNew('up', DG.TYPE.STRING);

  let addMilliseconds = function (time, ms) {
    return dayjs(time).add(ms, 'millisecond').toDate();
  };

  let tick = 20; // ms
  let start = addMilliseconds(new Date(), -tick * count);
  let next = function (symbol) {
    let v = tickers[symbol] + (Math.random() - 0.5) * 10;
    tickers[symbol] = v;
    return v;
  };
  let upDown = () => Math.random() > 0.5 ? '▲' : '▼';
  let getValues = (symbol, time) => [symbol, dayjs(time).toDate().toLocaleString('en-US'), next(symbol), Math.round(Math.random() * 10000), upDown()];

  let addTick = function (time) {
    for (let symbol in tickers)
      df.rows.addNew(getValues(symbol, time));
  };

  for (let i = 0; i < count; i++)
    addTick(addMilliseconds(start, tick));

  let view = grok.shell.addTableView(df);

  view.boxPlot({
    categoryColumnName: 'symbol',
    valueColumnName: 'price'
  });

  view.scatterPlot({
    xColumnName: 'time',
    yColumnName: 'price',
    sizeColumnName: 'volume',
    colorColumnName: 'symbol'
  });

  view.barChart({
    splitColumnName: 'symbol',
    valueColumnName: 'price',
    valueAggrType: 'avg'
  });

  view.grid.col('price').cellType = "percent completed";

  let row = 0;
  let timer = setInterval(function () {
    for (let symbol in tickers) {
      df.rows.setValues(Math.trunc(row % (count * Object.keys(tickers).length)),
        getValues(symbol, addMilliseconds(start, tick)));
      row++;
    }
    df.fireValuesChanged();
  }, tick);

  grok.events.onViewRemoved.subscribe(function (v) {
    if (v.name === view.name)
      clearInterval(timer);
  });
}

runStockBroker();
