//tags: DataQuery
//help-url: https://datagrok.ai/help/access/databases#parameterized-queries
// Manual master-details linking of tables that are dynamically retrieved from the database
grok.data.query('Samples:Countries', {}).then((countries) => {
  let customersView = null;
  grok.shell.addTableView(countries);
  countries.onCurrentRowChanged.subscribe((_) => {
    grok.data.query('Samples:CustomersInCountry', {country: countries.currentRow['country']}).then((t) => {
      //if (customersView === null)
        customersView = grok.shell.addTableView(t);
      customersView.dataFrame = t;
    });
    //grok.shell.addTableView(t);
  })
});

