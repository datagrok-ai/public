//api: DG.Stats.cumSum, DG.Stats.movingAvg
// Running totals and moving averages: values that depend on the rows before the current one.

let t = DG.DataFrame.fromCsv(
`region,date,amount
East,2024-01-03,10
West,2024-01-01,20
East,2024-01-02,
West,2024-01-04,40
East,2024-01-05,25`);
let amount = t.col('amount');

// In the natural row order. An empty value stays empty and does not break the total.
t.columns.add(amount.stats.cumSum()).name = 'running total';

// Along another column, restarting for each group.
t.columns.add(amount.stats.cumSum({orderBy: ['date'], by: ['region']})).name = 'by region, by date';

// Newest first.
t.columns.add(amount.stats.cumSum({orderBy: ['date'], ascending: [false]})).name = 'newest first';

// Trailing average over the last two rows by date. With minPeriods the first rows stay empty until the window is full.
t.columns.add(amount.stats.movingAvg(2, {orderBy: ['date']})).name = 'moving avg';
t.columns.add(amount.stats.movingAvg(2, {orderBy: ['date'], minPeriods: 2})).name = 'moving avg, full windows';

// Only the rows of a mask take part; the others stay empty.
let east = DG.BitSet.create(t.rowCount, (i) => t.get('region', i) === 'East');
t.columns.add(DG.Stats.fromColumn(amount, east).cumSum()).name = 'East only';

// An explicit visit order, such as the current sort order of a grid: grid.getRowOrder().
t.columns.add(amount.stats.cumSum({order: [4, 3, 2, 1, 0]})).name = 'bottom up';

grok.shell.addTableView(t);
