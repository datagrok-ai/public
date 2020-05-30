// Loading a csv file from the specified URL
grok.data.loadTable('https://public.datagrok.ai/demo/demog.csv')
  .then(t => grok.shell.addTableView(t));
