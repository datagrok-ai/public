// Providing custom comparer function

let months = { Jan: 1, Feb: 2, Mar: 3, Apr: 4, May: 5 };
let column = new DG.Column.fromList(DG.TYPE.STRING, 'months', ['Feb', 'Jan', 'May', 'Mar']);
column.valueComparer = (s1, s2) => months[s1] - months[s1];

let order = column.getSortedOrder();
grok.shell.info(order.map((i) => column.get(i)).join(', '));