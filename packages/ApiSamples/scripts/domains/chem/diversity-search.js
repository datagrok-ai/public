// An example of using diversity search.
// https://datagrok.ai/help/domains/chem/diversity-search

let mols = await grok.chem.diversitySearch(grok.data.testData('molecules', 100).col('smiles'));
grok.shell.addTableView(mols);