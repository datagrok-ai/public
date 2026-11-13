//api: DG.ViewLayout.boundColumnNames
// The columns a stored layout's viewers are bound to, the Filters viewer's filters included,
// read from the layout itself: no view has to show it, and a view state works as well as a saved layout.

const tv = grok.shell.addTableView(grok.data.demo.demog());
tv.scatterPlot({x: 'height', y: 'weight', color: 'race'});
tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'sex'});

const layout = DG.ViewLayout.fromViewState(tv.saveLayout().viewState);
grok.shell.info(`Bound columns: ${layout.boundColumnNames.join(', ')}`);
