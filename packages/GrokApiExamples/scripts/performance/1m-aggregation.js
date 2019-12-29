// Aggregate a table with 1 million rows
//
// See also: https://public.datagrok.ai/help/concepts/performance

let wells = time('create', () => grok.testData('wells', 1000000));

let concentrations = time('aggregate', () => wells
    .groupBy(['row', 'role'])
    .avg('concentration')
    .aggregate());

grok.addTableView(concentrations);