grok.functions
  .eval('OpenServerFile("System:DemoFiles/beer.csv")')
  .then((t) => grok.shell.addTableView(t[0]));
