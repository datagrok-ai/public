// An example of using substructure search.

grok.data.loadDataFrame('https://public.datagrok.ai/demo/sar_small.csv')
    .then(t => grok.chem.substructureSearch(t.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1', false)
        .then(function (bs) {
            t.selection.copyFrom(bs);
            grok.shell.addTableView(t);
        }));
