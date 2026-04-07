/**
 * Opening a table with one million columns.
 * See also: https://datagrok.ai/help/develop/performance
 */

let table = grok.data.testData('random walk', 10, 1000000);

grok.shell.addTableView(table);
