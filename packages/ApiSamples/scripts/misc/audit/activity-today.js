const view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:TopFunctions', {date: 'today', users: ['all']})
  .then((t) => view.append(ui.block(DG.Viewer.scatterPlot(t).root)));