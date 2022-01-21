//name: AddMorganFPCache
//description: Adds a column with Morgan fingerprints
//language: javascript

(async () => {
  let inDf = grok.shell.tables[0];
  let molCol = inDf.columns[0];
  const length = molCol.length;
  let fpCol = await chem.getMorganFingerprints(molCol);
  let morganFpCol = DG.Column.fromType(DG.TYPE.STRING, 'cache_morganFP', length);
  for (let i = 0; i < length; ++i) {
    morganFpCol.set(i, fpCol.get(i).toBinaryString());
  }
  inDf.columns.add(morganFpCol);
})();