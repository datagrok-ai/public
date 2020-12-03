//name: testSubstructureSearch
//description: Substructure search via RDKit
//language: javascript

async function time(name, n, f) {
  let start = new Date();
  for (let k = 0; k < n; ++k) {
    await f();
  }
  let stop = new Date();
  console.log(`${name}: ${(stop - start) / n} ms`);
}

(async () => {
  
  let searchFor = ['c1ccccc1', 'C1CCCCC1', 'CC'];
  
  let df = await grok.data.getDemoTable('chem/zbb/99_p3_4.5-6.csv');
  const N = Math.min(2000, df.rowCount);
  df.rows.removeAt(N, df.rowCount - N, false);
  
  const n = 10;
  
  await time('Substructure Search - Graph', n, async () => {
    for (let s of searchFor) {
      let result = await grok.chem.substructureSearch(df.col('smiles'), s, { substructLibrary: false });
      // console.log(s + ": " + result.toString());
    }
  });
  
  await time('Substructure Search - Library - Init', 1, async() => {
    await grok.chem.substructureSearch(df.col('smiles'));
  });
  
  await time('Substructure Search - Library - Search', n, async () => {
    for (let s of searchFor) {
      let result = await grok.chem.substructureSearch(df.col('smiles'), s);
      // console.log(s + ": " + result.toString());
    }
  });
  
})();