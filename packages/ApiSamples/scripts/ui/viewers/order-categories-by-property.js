// Order a categorical axis by a number computed from the values themselves.
// A `categoryOrderer` function receives the distinct values of its semantic type and
// returns one comparable key per value; viewers offer it in the axis menu's Order... dialog.
// An optional second string parameter with `choices` turns one function into one property
// per choice, as in `Chem:molecularProperty(MW)`.
const f = grok.functions.register({
  signature: 'column nameLengthKeys(column/DemoNames values)',
  tags: 'categoryOrderer',
  run: (values) => {
    const col = DG.toJs(values);
    return DG.Column.fromList('double', 'keys', col.toList().map((s) => (s ?? '').length));
  },
});

const df = DG.DataFrame.fromCsv('name,value\nviolet,1\nred,2\ncyan,3\nred,4\nviolet,5');
df.col('name').semType = 'DemoNames';
const view = grok.shell.addTableView(df);

// the order is a JSON string with the same id a user would pick from the menu;
// also: yCategoryOrder (scatterplot), xCategoryOrder (line chart), categoryOrder (box plot; bar and pie charts
// under the 'custom' sort type), xCategoryOrder / yCategoryOrder / zCategoryOrder (3D scatter plot)
view.scatterPlot({x: 'name', y: 'value'}).setOptions({xCategoryOrder: JSON.stringify({property: f.nqName, order: 'asc'})});
