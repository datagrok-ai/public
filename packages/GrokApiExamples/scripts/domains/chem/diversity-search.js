// An example of using diversity search.
//
// https://datagrok.ai/help/domains/chem/diversity-search

grok.chem.diversitySearch(grok.testData('molecules', 100).col('smiles'))
    .then(function (mols) {
        let col = grok.Column.fromStrings('smiles', mols);
        grok.addTableView(DataFrame.fromColumns([col]));
    });
