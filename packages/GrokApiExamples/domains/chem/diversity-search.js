// An example of using diversity search.
//
// https://datagrok.ai/help/domains/chem/diversity-search

chem.diversitySearch(gr.testData('molecules', 100).col('smiles'))
    .then(function (mols) {
        let col = Column.fromStrings('smiles', mols);
        gr.addTableView(DataFrame.fromColumns([col]));
    });
