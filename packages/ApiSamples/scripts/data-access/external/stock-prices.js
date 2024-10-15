grok.data.loadTable('https://calmcode.io/static/data/bigmac.csv')
  .then(t => grok.shell.addTableView(t));
