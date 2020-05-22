// An example of using diversity search.
//
// https://datagrok.ai/help/domains/chem/diversity-search

grok.chem.diversitySearch(grok.data.testData('molecules', 100).col('smiles'))
    .then(function (mols) {
        let col = DG.Column.fromStrings('smiles', mols);
        grok.shell.addTableView(DG.DataFrame.fromColumns([col]));
    });
