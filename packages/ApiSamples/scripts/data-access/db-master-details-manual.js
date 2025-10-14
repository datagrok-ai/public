//tags: DataQuery
//help-url: https://datagrok.ai/help/access/databases#parameterized-queries
// Manual master-details linking of tables that are dynamically retrieved from the database
let countries = await grok.data.query('Samples:Countries', {})
let customersView = null;

grok.shell.addTableView(countries);
countries.onCurrentRowChanged.subscribe(async (_) => {
  let t = await grok.data.query('Samples:CustomersInCountry', {country: countries.currentRow['country']})
  customersView = grok.shell.addTableView(t);
  customersView.dataFrame = t;
})

