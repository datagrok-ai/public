//tags: ViewLayout, Dataframe
let df = grok.data.demo.demog();
let view = grok.shell.addTableView(df);
grok.dapi.layouts.getApplicable(df).then(layouts => {
    if (layouts.length > 0)
        view.loadLayout(layouts[0]);
});
