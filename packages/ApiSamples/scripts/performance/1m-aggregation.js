/**
 * Aggregate a table with 1 million rows
 * See also: https://datagrok.ai/help/develop/performance
 */

let time = function(s, f) {
  let start = new Date();
  let result = f();
  let stop = new Date();
  console.log(`${s}: ${stop - start}ms`);
  grok.shell.info(`${s}: ${stop - start}ms`);
  return result;
};

let wells = time('create', () => grok.data.testData('wells', 1000000));

let concentrations = time('aggregate', () => wells
  .groupBy(['row', 'role'])
  .avg('concentration')
  .aggregate());

grok.shell.addTableView(concentrations);
