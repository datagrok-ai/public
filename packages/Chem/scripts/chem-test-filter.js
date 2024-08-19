//name: chemTestFilter
//description: Test a packaged molecular filter on a selected dataframe
//language: javascript

(async () => {
  const df = grok.data.demo.molecules(1000);
  await grok.data.detectSemanticTypes(df);
  const filter = await grok.functions.call('Chem:substructureFilter');
  filter.attach(df);
  grok.shell.addTableView(df);
  const colChoice = ui.input.column('Column', {value: filter.column, table: filter.dataFrame, filter: (col) => {
    filter.column = col;
    filter.dataFrame.filter.setAll(true, false);
    filter.dataFrame.rows.requestFilter();
  }});
  ui.dialog({title: 'Chem Filter'})
    .add(colChoice)
    .add(filter.root)
    .show();
})();
