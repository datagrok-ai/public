// Opening a table with one hundred million rows.
//
// We can do one billion rows as well, but it will take about two minutes to execute.
// Note that viewers that visualize each point (such as scatter plot) might become unresponsive,
// but will still render the picture after all.
//
// See also: https://public.datagrok.ai/help/concepts/performance

let table = gr.testData('random walk', 100000000, 2);

gr.addTableView(table);