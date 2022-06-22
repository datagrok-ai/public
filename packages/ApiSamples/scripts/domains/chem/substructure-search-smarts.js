//name: testSubstructureSearchSmarts
//language: javascript

// An example of using substructure search with SMARTS pattern

grok.data.getDemoTable('chem/smiles.csv')
  .then(t => grok.chem.searchSubstructure(
    t.col('canonical_smiles'), '[!#6&!#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1')
      .then(function (bs) {
        t.filter.copyFrom(bs);
        grok.shell.addTableView(t);
      }));