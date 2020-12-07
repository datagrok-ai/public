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
  const N = Math.min(20000, df.rowCount);
  df.rows.removeAt(N, df.rowCount - N, false);
  
  const q = 3;
  
  console.log('Similarity scoring benchmark');
  
  await time(`Building library for ${df.rowCount} molecules`, 1, async() => {    
    await grok.chem.similarityScoring(
      df.col('smiles'), df.col('smiles').get(0));
  });
  
  await time(`Searching the first ${q} molecules`, 3, async() => {    
    for (let i = 0; i < q; ++i) {
      let similar = await grok.chem.similarityScoring(
        df.col('smiles'), df.col('smiles').get(i));
    }
  });
  
})();