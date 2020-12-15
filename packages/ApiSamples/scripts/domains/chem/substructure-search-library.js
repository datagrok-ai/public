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
  df.rows.removeAt(N - 1, df.rowCount - N + 1, false);
  let col = df.col('smiles');
  
  const n = 5;
  
  console.log('Substructure search microbenchmark');
  
  await time(`Graph-based search on ${df.rowCount} molecules`, n, async () => {
    for (let s of searchFor) {
      let result = await grok.chem.substructureSearch(col, s, { substructLibrary: false });
      // console.log(s + ": " + result.toString());
    }
  });
  
  await time(`Building a library for ${df.rowCount} molecules`, 1, async() => {
    await grok.chem.substructureSearch(col, '');
  });
  
  await time(`Searching the library on ${searchFor.length} molecules`, n, async () => {
    for (let s of searchFor) {
      let result = await grok.chem.substructureSearch(col, s);
      // console.log(s + ": " + result.toString());
    }
  });
  
})();