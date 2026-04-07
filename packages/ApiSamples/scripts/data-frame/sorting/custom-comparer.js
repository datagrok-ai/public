// Providing custom comparer function

// Let's say we have a string column where values have some logical order
let column = DG.Column.fromList(DG.TYPE.STRING, 'months', ['Feb', 'Jan', 'May', 'Mar']);

// Define and set a custom value comparer that is used for sorting later
let months = {Jan: 1, Feb: 2, Mar: 3, Apr: 4, May: 5};
column.valueComparer = (s1, s2) => months[s1] - months[s2];

// It is possible to get the sorted order if needed
let order = column.getSortedOrder();
grok.shell.info(Array.from(order).map((i) => column.get(i)).join(', '));

// Or, simply double-click on the header to sort it visually
grok.shell.addTableView(DG.DataFrame.fromColumns([column]));