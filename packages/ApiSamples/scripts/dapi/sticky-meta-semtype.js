//name: sticky-meta-semtype
//tags: demo
//language: javascript

// uniqueName ensures that sample creates new configuration
var uniqueName = (prefix) => 'apisamles-' + prefix + '-' + (Math.random() + 1).toString(36).substring(7);

// create sample dataset on molecules
var df = grok.data.testData('molecules', 15);

// Create a schema named Date, connected with molecules, and having 1 parameter "date"
var schema = await grok.dapi.stickyMeta.createSchema(
  uniqueName('Date'),
  [{name: uniqueName('Molecule'), matchBy: 'semtype=molecule'}],
  [{name: 'date', type: 'datetime'}]
);

// Open dataframe in a view.
var startView = grok.shell.addTableView(df);

// This timeout is needed for Chem package to apply its semtype detector on opened view
setTimeout(async () => {

  // Create sample date time values column with one value
  var valuesColumn = DG.Column.dateTime('date', 15);
  valuesColumn.set(0, '03/20/2024');
  
  // Save date as sticky meta. 
  await grok.dapi.stickyMeta.setAllValues(schema, df.columns.byName('smiles'), DG.DataFrame.fromColumns([valuesColumn]));
  
  // Trigger cell redraw to fetch latest changes
  // This line can be removed and replaced by manual scroll of grid or column resize
  startView.grid.invalidate();
}, 1000);

