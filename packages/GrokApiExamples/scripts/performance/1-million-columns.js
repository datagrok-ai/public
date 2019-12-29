// Opening a table with one million columns.
//
// See also: https://public.datagrok.ai/help/concepts/performance

let table = grok.testData('random walk', 10, 1000000);

grok.addTableView(table);