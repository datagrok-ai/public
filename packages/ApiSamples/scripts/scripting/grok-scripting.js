grok.functions
  .eval('OpenServerFile("Skalkin:DatagrokData/formats/xslx/beer.xlsx")')
  .then((t) => grok.shell.addTableView(DG.toJs(t[0])));