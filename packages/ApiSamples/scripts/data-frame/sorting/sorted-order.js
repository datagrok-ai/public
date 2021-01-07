// DataFrame.getSortedOrder does not sort rows, but instead returns the sorted order,
// which could be used further for multiple purposes.

let t = grok.data.demo.demog();
let order = t.getSortedOrder(['race', 'age']);
grok.shell.info(`min: ${t.columns.byName('race').get(order[0])}`);
grok.shell.info(`min: ${t.columns.byName('race').get(order[order.length - 1])}`);