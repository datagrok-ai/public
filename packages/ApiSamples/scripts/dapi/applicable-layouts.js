//tags: ViewLayout, Dataframe
let df = grok.data.demo.demog();
let view = grok.shell.addTableView(df);
let layouts = await grok.dapi.layouts.getApplicable(df);
if (layouts.length > 0) 
    view.loadLayout(layouts[0]);
