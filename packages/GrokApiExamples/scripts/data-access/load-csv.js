// Loading a csv file from the specified URL

grok.loadDataFrame('/demo/demog.csv')
  .then(t => grok.addTableView(t));
