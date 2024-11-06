//help-url: https://datagrok.ai/help/govern/user
const view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:UniqueUsersList', {date: 'yesterday'})
  .then((t) => {
    const ids = new Set(t.getCol('id').toList());
    grok.dapi.getEntities([...ids]).then((users) => view.append(ui.list(users)));
  });