// adding the specified filter to the already existing filter group

let tv = grok.shell.addTableView(grok.data.demo.demog(100));
let options = {
  type: DG.FILTER_TYPE.CATEGORICAL,
  column: 'site',
  selected: ['New York', 'Orlando']
};
tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd(options);