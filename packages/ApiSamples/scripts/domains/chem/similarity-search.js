// An example of using similarity search.
//
// https://datagrok.ai/help/domains/chem/similarity-search

grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv')
    .then(molecules => grok.chem.similaritySearch(molecules.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1')
        .then(similar => grok.shell.addTableView(similar)));

