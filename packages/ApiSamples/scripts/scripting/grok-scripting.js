let t = await grok.functions.eval('OpenServerFile("System:DemoFiles/beer.csv")');
grok.shell.addTableView(t[0]);
