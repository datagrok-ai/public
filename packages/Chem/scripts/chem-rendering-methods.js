//name: chemBenchmark
//language: javascript

function randomInt(min, max) {
  return min + Math.floor((max - min) * Math.random());
}

/* TODO: internalize with Datagrok */
async function time(name, n, f) {
  let start = window.performance.now();
  for (let k = 0; k < n; ++k) {
    await f();
  }
  let stop = window.performance.now();
  console.log(`${name}: ${((stop - start) / n).toFixed(2)} ms`);
}

async function timeBySubarrays(name, samples, sz, n, f) {
  console.assert(sz <= samples.length,
    'Subsample size is larger than samples array');
  let total = 0;
  let cntr = n;
  while (cntr > 0) {
    let from = randomInt(0, samples.length - 1 - sz);
    let to = from + sz;
    let sample = samples.slice(from, to);
    const start = window.performance.now();
    await f(sample);
    const stop = window.performance.now();
    total += stop - start;
    cntr--;
  }
  console.log(
    `${name}, averaged on ${n} subarrays of length ` +
    `${sz} from a pool of ${samples.length} items: ${total.toFixed(2)} ms`);
}

function getRandomSubarray(arr, size) {
  var shuffled = arr.slice(0), i = arr.length, temp, index;
  while (i--) {
    index = Math.floor((i + 1) * Math.random());
    temp = shuffled[index];
    shuffled[index] = shuffled[i];
    shuffled[i] = temp;
  }
  return shuffled.slice(0, size);
}

/* TODO: make this uniform with timeBySubarrays */
async function timeByRandomArrays(name, samples, sz, n, f) {
  console.assert(sz <= samples.length,
    'Subsample size is larger than samples array');
  let total = 0;
  let cntr = n;
  while (cntr > 0) {
    let sample = getRandomSubarray(samples, sz);
    const start = window.performance.now();
    await f(sample);
    const stop = window.performance.now();
    total += stop - start;
    cntr--;
  }
  console.log(
    `${name}, averaged on ${n} random subsets of size ` +
    `${sz} from a pool of ${samples.length} items: ${total.toFixed(2)} ms`);
}

function createCanvas(width, height) {
  let c = document.createElement('canvas');
  c.setAttribute('width', width);
  c.setAttribute('height', height);
  return c;
}

(async () => {

  const N = 1000;
  const n = 20;
  const t = 100;
  const x = 0;
  const y = 0;
  const w = 200;
  const h = 100;
  const mu = 5;

  let df = await grok.data.getDemoTable('chem/zbb/99_p3_4.5-6.csv');
  let col = df.col('smiles');

  console.log('Chem Benchmark');

  let canvas = createCanvas(w, h);
  let module = _ChemPackage.rdKitModule;
  let molStrings = col.toList().slice(0, N - 1); // TODO: add and try random with a seed
  let molArray = molStrings.map(s => module.get_mol(s));

  const renderMol = (mol) => {
    const opts = {
      "clearBackground": false,
      "offsetx": Math.floor(x),
      "offsety": -Math.floor(y),
      "width": Math.floor(w),
      "height": Math.floor(h)
    };
    mol.draw_to_canvas_with_highlights(
      canvas, JSON.stringify(opts));
  };

  /* --- */
  console.log('1. Rendering pre-built molecules...'); /* --- */

  const renderMols = (mols) => {
    for (let mol of mols) {
      renderMol(mol);
    }
  };

  await time(`Rendering a ${N} molecules, pre-built`, 1, async () => {
    renderMols(molArray);
  });

  /* --- */
  console.log('2. Rendering molecules, without LRU-cache...'); /* --- */

  const renderMolsByStrings = (strings) => {
    for (let smiles of strings) {
      let mol = module.get_mol(smiles);
      renderMol(mol);
      mol.delete();
    }
  };

  await time(`Rendering a ${N} molecules, no cache`, 1, async () => {
    renderMolsByStrings(molStrings);
  });

  /* --- */
  console.log('3. Rendering molecules, with LRU-cache...'); /* --- */

  let rendererCache = new DG.LruCache();
  rendererCache.onItemEvicted = (mol) => mol.delete();
  const renderMolsWithCache = (strings) => {
    for (let smiles of strings) {
      let mol = rendererCache.getOrCreate(smiles, s => module.get_mol(s));
      renderMol(mol);
    }
  };

  await time(`Rendering a ${N} molecules, with cache`, 1, async () => {
    renderMolsWithCache(molStrings);
  });

  /* --- */
  console.log('4. Horizontal scrolling...'); /* --- */

  await timeBySubarrays(`Rendering ${n} molecules ${t} times`, molStrings, n, 5, (sample) => {
    for (let j = 0; j < t; ++j) {
      renderMolsWithCache(sample);
    }
  });

  /* --- */
  console.log('5. Vertical scrolling...'); /* --- */

  await timeBySubarrays(`Rendering ${n} molecules ${t} times with a sliding window`, molStrings, n + t, 5, (sample) => {
    for (let j = 0; j < t; ++j) {
      renderMolsWithCache(sample.slice(j, j + n));
    }
  });

})();