let tv = grok.shell.addTableView(grok.data.demo.demog());
tv.filters({filters: [
    {type: 'Widgets:radioButtonFilter', columnName: 'race'},
  ]});