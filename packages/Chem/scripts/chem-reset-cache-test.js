//name: resetCacheTest
//language: javascript

(async () => {

  const sample = 'COc1ccc2cc(ccc2c1)C(C)C(=O)N3CCCC3C(=O)OCCCc4cccnc4';
  let df = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'smiles', [
      'CC(C)Cc1ccc(cc1)C(C)C(=O)OCCCc2cccnc2',
      'CC(C(=O)OCCCc1cccnc1)c2cccc(c2)C(=O)c3ccccc3',
      'COc1ccc(NC(=O)NN2C(=Nc3ccccc3C2=O)C)cc1'
    ])
  ]);
  let col = df.col('smiles');
  let res1 = await grok.chem.getSimilarities(col);
  let resCol1 = await grok.chem.getSimilarities(col, sample);
  grok.shell.addTableView(DG.DataFrame.fromColumns([resCol1]));
  df.rows.removeAt(0, 2);
  df.rows.addNew(['Cc1noc(n1)c2ccccn2']);
  df.rows.addNew(['Cc1c(COc2cccc(F)c2)oc3cccc(OCCNCc4cccnc4)c13']);
  df.rows.addNew(['CCCCCCCCCCCCCCCC(=O)OCCCCCCOC(=O)CCCCCCCCCCCCCCC']);
  let resCol2 = await grok.chem.getSimilarities(col, sample);
  grok.shell.addTableView(DG.DataFrame.fromColumns([resCol2]));

})();