// An example of using diversity search.
//
// https://datagrok.ai/help/domains/chem/diversity-search

chem.diversitySearch(grok.testData('molecules', 100).col('smiles'))
    .then(function (mols) {
        let col = Column.fromStrings('smiles', mols);
        grok.addTableView(DataFrame.fromColumns([col]));
    });
