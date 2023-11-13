grok.functions
  .eval('OpenServerFile("Samples:Files/beer.csv")')
  .then((t) => grok.shell.addTableView(t[0]));
