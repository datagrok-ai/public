//help-url: https://datagrok.ai/help/govern/user
const view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:Usage', {date: 'today', users: ['all']})
  .then((t) => {
    const ids = new Set(t.getCol('user_id').toList());
    grok.dapi.getEntities([...ids]).then((users) => view.append(ui.list(users)));
  });