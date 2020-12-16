//name: chemBenchmark
//language: javascript

// TODO: internalize with Datagrok
async function time(name, n, f) {
  let start = new Date();
  for (let k = 0; k < n; ++k) {
    await f();
  }
  let stop = new Date();
  console.log(`${name}: ${(stop - start) / n} ms`);
}

function createCanvas(width, height) {
  var c = document.createElement('canvas');
  c.setAttribute('width', width);
  c.setAttribute('height', height);
  return c;
}

(async () => {

  const N = 1000;
  const x = 0;
  const y = 0;
  const w = 200;
  const h = 100;
  
  let df = await grok.data.getDemoTable('chem/zbb/99_p3_4.5-6.csv');
  if (N < df.rowCount)
  	df.rows.removeAt(N, df.rowCount - N, false);
  let col = df.col('smiles');
  
  console.log('Chem Benchmark');
  
  let canvas = createCanvas(w, h);
  let module = _ChemPackage.rdKitModule;
  let molStrings = col.toList();
  let molArray = molStrings.map(s => module.get_mol(s));
  
  await time(`Rendering a ${N} molecules`, 1, async () => {
    for (let mol of molArray) {
   	  const opts = {
      	"clearBackground": false,
        "offsetx": Math.floor(x),
        "offsety": -Math.floor(y),
        "width": Math.floor(w),
        "height": Math.floor(h)
      };
      mol.draw_to_canvas_with_highlights(
        canvas, JSON.stringify(opts));
    }
  });
  
  molArray.forEach(m => m.delete());
  
})();