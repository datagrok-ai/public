let tv = grok.shell.addTableView(grok.data.demo.demog());

// We need to load the external package before we use the
// synchronous "tv.filters" method that uses the filter
grok.functions.call('Widgets:radioButtonFilter').then((_) => {
  tv.filters({filters: [
      {type: 'Widgets:radioButtonFilter', columnName: 'race'},
    ]});
});