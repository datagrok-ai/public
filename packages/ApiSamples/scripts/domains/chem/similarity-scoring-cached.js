//name: testSimilarityScoringCached
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
  
  let df = await grok.data.getDemoTable('chem/zbb/99_p3_4.5-6.csv');
  const N = Math.min(2000, df.rowCount);
  df.rows.removeAt(N - 1, df.rowCount - N + 1, false);
  let col = df.col('smiles');
  
  const q = 3;
  const n = 5;
  
  console.log('Similarity scoring microbenchmark');
  
  await time(`Building library for ${df.rowCount} molecules`, 1, async() => {    
    await grok.chem.similarityScoring(col, '');
  });
  
  await time(`Searching the first ${q} molecules`, n, async() => {    
    for (let i = 0; i < q; ++i) {
      let similar = await grok.chem.similarityScoring(col, col.get(i));
    }
  });
  
})();