let table = DG.DataFrame.fromCsv(
`country, jan, feb, mar, apr, may
USA, 23, 34, 56, 5, 7
Canada, 2, 3, 5, 4, 5`);

let unpivoted = table.unpivot(['country'], ['jan', 'feb', 'mar', 'apr', 'may']);

grok.shell.add(unpivoted);