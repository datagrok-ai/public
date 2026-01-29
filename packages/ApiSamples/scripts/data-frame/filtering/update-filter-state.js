//updating state of opened filter
let beer = await grok.data.getDemoTable('beer.csv');
let view = grok.shell.addTableView(beer);

setTimeout(function() {
  let filters = view.filters();
  filters.updateOrAdd({
    type: 'text',
    column: 'Flavor',
    gridNames: ['malt'],
    value: 'sour'
  });
}, 500);