//tags: ViewLayout, Dataframe
let df = grok.data.demo.demog();
let view = grok.shell.addTableView(df);
grok.dapi.layouts.getApplicable(df).then(layouts => view.loadLayout(layouts[0]));
