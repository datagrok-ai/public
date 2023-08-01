grok.functions
  .eval('OpenServerFile("Demo:Files/beer.csv")')
  .then((t) => grok.shell.addTableView(t[0]));
