//name: sticky-meta-2-semtype
//language: javascript

// uniqueName ensures that created sticky meta configuration will be unique.
let uniqueName = (prefix) => 'apisamples-' + prefix + '-' + (Math.random() + 1).toString(36).substring(7);

// Create sample dataset of molecules.
let df = grok.data.testData('molecules', 3);

// Create a schema named Date, connected with molecules, and having 1 parameter "int-meta".
let schema = await grok.dapi.stickyMeta.createSchema(
  uniqueName('int-meta'),
  [{name: uniqueName('Molecule'), matchBy: 'semtype=molecule'}],
  [{name: 'int-meta', type: 'int'}]
);

// Open dataframe in a view.
let startView = grok.shell.addTableView(df);

// We wait for Chem package to apply its semtype detector on an opened dataframe.
df.onSemanticTypeDetected.subscribe(async (_) => {

  // Create sample values column. Name of the column matches name of the property.
  let valuesColumn = DG.Column.fromList('int', 'int-meta', [1, 2, 3]);

  // Save date as sticky meta. 
  await grok.dapi.stickyMeta.setAllValues(schema, df.columns.byName('smiles'), DG.DataFrame.fromColumns([valuesColumn]));

  // Trigger cell redraw to redraw grid so grid shows correct sticky meta indication.
  // This action can be replaced by manual scroll of grid or column width resize.
  startView.grid.invalidate();
});
