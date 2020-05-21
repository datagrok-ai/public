// The following demo shows, how to work with dynamically changing data in Datagrok.

var count = 50;
var tickers = new Map([
    ['AAPL', 200],
    ['TSLA', 200],
    ['GOOG', 200],
    ['NFLX', 200]
]);

var df = grok.DataFrame.create(0);
df.name = 'stocks';
df.columns.addNew('symbol', grok.TYPE_STRING);
df.columns.addNew('time', grok.TYPE_DATE_TIME);
df.columns.addNew('price', grok.TYPE_FLOAT);
df.columns.addNew('volume', grok.TYPE_FLOAT);
df.columns.addNew('up', grok.TYPE_STRING);

var addMilliseconds = function (time, ms) {
    time.setMilliseconds(time.getMilliseconds() + ms);
    return time;
};

var tick = 20; // ms
var start = addMilliseconds(new Date(), -tick * count);
let next = function (symbol) {
    let v = tickers.get(symbol) + (Math.random() - 0.5) * 10;
    tickers.set(symbol, v);
    return v;
};
var upDown = () => Math.random() > 0.5 ? '▲' : '▼';
var getValues = (symbol, time) => [symbol, DateTime.fromDate(time).d, next(symbol), Math.round(Math.random() * 10000), upDown()];

let addTick = function (time) {
    for (let symbol of tickers.keys())
        df.rows.addNew(getValues(symbol, time));
};

for (let i = 0; i < count; i++)
    addTick(addMilliseconds(start, tick));

var view = grok.addTableView(df);

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

var row = 0;
var timer = setInterval(function() {
    for (let symbol of tickers.keys()) {
        df.rows.setValues(Math.trunc(row % (count * tickers.size)), getValues(symbol, addMilliseconds(start, tick)));
        row++;
    }
    df.fireValuesChanged();
}, tick);

grok.onViewRemoved(function (v) {
    if (v.name === view.name)
        clearInterval(timer);
});
