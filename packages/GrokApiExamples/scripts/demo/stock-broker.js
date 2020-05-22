// The following demo shows, how to work with dynamically changing data in Datagrok.

let count = 50;
let tickers = new Map([
    ['AAPL', 200],
    ['TSLA', 200],
    ['GOOG', 200],
    ['NFLX', 200]
]);

let df = DG.DataFrame.create(0);
df.name = 'stocks';
df.columns.addNew('symbol', grok.TYPE_STRING);
df.columns.addNew('time', grok.TYPE_DATE_TIME);
df.columns.addNew('price', grok.TYPE_FLOAT);
df.columns.addNew('volume', grok.TYPE_FLOAT);
df.columns.addNew('up', grok.TYPE_STRING);

let addMilliseconds = function (time, ms) {
    time.setMilliseconds(time.getMilliseconds() + ms);
    return time;
};

let tick = 20; // ms
let start = addMilliseconds(new Date(), -tick * count);
let next = function (symbol) {
    let v = tickers.get(symbol) + (Math.random() - 0.5) * 10;
    tickers.set(symbol, v);
    return v;
};
let upDown = () => Math.random() > 0.5 ? '▲' : '▼';
let getValues = (symbol, time) => [symbol, DateTime.fromDate(time).d, next(symbol), Math.round(Math.random() * 10000), upDown()];

let addTick = function (time) {
    for (let symbol of tickers.keys())
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
let timer = setInterval(function() {
    for (let symbol of tickers.keys()) {
        df.rows.setValues(Math.trunc(row % (count * tickers.size)), getValues(symbol, addMilliseconds(start, tick)));
        row++;
    }
    df.fireValuesChanged();
}, tick);

grok.events.onViewRemoved.subscribe(function (v) {
    if (v.name === view.name)
        clearInterval(timer);
});
