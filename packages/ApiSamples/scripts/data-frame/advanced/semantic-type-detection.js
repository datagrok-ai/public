// Custom handling of the semantic type detection
// See https://datagrok.ai/help/discover/semantic-types

let table = grok.data.parseCsv('smiles\nFc1cc(Cl)ccc1Br');

// Normally, the smiles column will be detected as 'Molecule' when we create a table view.
// However, this time we will override it

table.onSemanticTypeDetecting.subscribe((_) => {
  table.col('smiles').semType = 'bananas';
});

table.onSemanticTypeDetected.subscribe((_) => {
  grok.shell.info(table.col('smiles').semType);
});

grok.shell.addTableView(table);
