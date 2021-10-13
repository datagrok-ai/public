//name: chemTestFilter
//description: Test a packaged molecular filter on a selected dataframe
//language: javascript

grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv').then(
  async (df) => {
    await grok.data.detectSemanticTypes(df);
    let filter = await grok.functions.call("Chem:substructureFilter");
    filter.attach(df);
    grok.shell.addTableView(df);
    let colChoice = ui.columnInput('Column', filter.dataFrame, filter.column, (col) => {
      filter.column = col;
      filter.dataFrame.filter.setAll(true, false);
      filter.dataFrame.rows.requestFilter();
    });
    ui.dialog({title: 'Chem Filter'})
      .add(colChoice)
      .add(filter.root)
      .show();
  }
)