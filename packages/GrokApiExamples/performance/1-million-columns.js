// Opening a table with one million columns.
//
// See also: https://public.datagrok.ai/help/concepts/performance

let table = gr.testData('random walk', 10, 1000000);

gr.addTableView(table);