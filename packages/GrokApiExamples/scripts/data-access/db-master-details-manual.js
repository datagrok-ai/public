// Manual master-details linking of tables that are dynamically retrieved from the database

grok.query('northwind:countries', {}).then((countries) => {
    var customersView = null;
    grok.addTableView(countries);
    countries.onCurrentRowChanged((_) => {
        grok.query('northwind:customersByCountry', {country: countries.currentRow['country']}).then((t) => {
            if (customersView == null)
                customersView = grok.addTableView(t);
            customersView.dataFrame = t;
        });
        grok.addTableView(countries);
    })
});