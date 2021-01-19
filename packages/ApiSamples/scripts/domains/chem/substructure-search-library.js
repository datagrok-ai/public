//name: testSubstructureSearch
//language: javascript

// An example of using substructure search.

grok.data.loadTable('https://public.datagrok.ai/demo/sar_small.csv')
  .then(t => grok.chem.searchSubstructure(t.col('smiles'), 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1')
    .then(function (bs) {
      t.filter.copyFrom(bs);
      grok.shell.addTableView(t);
    }));