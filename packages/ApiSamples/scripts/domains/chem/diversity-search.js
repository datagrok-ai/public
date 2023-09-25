// An example of using diversity search.
// https://datagrok.ai/help/domains/chem/diversity-search

grok.chem.diversitySearch(grok.data.testData('molecules', 100).col('smiles'))
  .then(function(mols) {
    grok.shell.addTableView(mols);
  });