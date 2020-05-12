// Loading a csv file from the specified URL

grok.loadDataFrame('https://public.datagrok.ai/demo/demog.csv')
  .then(t => grok.addTableView(t));
