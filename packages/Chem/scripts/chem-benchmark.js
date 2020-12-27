//name: chemBenchmark
//language: javascript

const randomInt = (min, max) => min + Math.floor((max - min) * Math.random());
const range = (l, r) => new Array(r - l).fill().map((_, k) => k + l);

// TODO: internalize with Datagrok
async function time(name, n, f) {
  let start = window.performance.now();
  for (let k = 0; k < n; ++k) {
    await f();
  }
  let stop = window.performance.now();
  console.log(`${name}: ${((stop - start) / n).toFixed(2)} ms`);
}

function getIdxRandomSubset(N, n) {
  let arr = range(0, N);
  let shuffled = arr.slice(0);
  let i = arr.length, temp, index;
  while (i--) {
    index = randomInt(0, i + 1);
    temp = shuffled[index];
    shuffled[index] = shuffled[i];
    shuffled[i] = temp;
  }
  return shuffled.slice(0, n);
}

function getIdxRandomSubarray(N, n) {
  const from = randomInt(0, N - 1 - n);
  const to = from + n;
  return range(from, to);
}

let caseCntr = 0;

async function testCase(title, f) {
  await time(++caseCntr + '. ' + title, 1, f);
}

(async () => {

  const nRender = 1000;
  const nScroll = 20;
  const tScroll = 100;
  const nSearch = 25000;
  const nSample = 10;
  const x = 0, y = 0;
  const w = 200, h = 100;

  let df = await grok.data.getDemoTable('chem/zbb/99_p3_4.5-6.csv');
  // For using with your HOME files in Datagrok:
  // let df = (await grok.functions.eval('OpenServerFile("UserName:Home/Chembl_100K.csv")'))[0];
  const colName = 'smiles';
  let col = df.col(colName);
  let N = col.length;

  console.log('Chem Benchmark');

  let rdKitCellRenderer = new RDKitCellRenderer();
  let grid = DG.Viewer.grid(df);
  let canvas = grid.canvas;
  let ctx = canvas.getContext('2d');
  let rowsInd = null;

  rowsInd = getIdxRandomSubset(N, nRender);
  await testCase(`Rendering a ${rowsInd.length} random molecules`, () => {
    for (i of rowsInd) {
      rdKitCellRenderer.render(ctx, x, y, w, h, grid.cell(colName, i), null);
    }
  });

  rowsInd = getIdxRandomSubarray(N, nScroll);
  await testCase(`Horizontal scrolling (${rowsInd.length} random molecules, ${tScroll} times)`, () => {
    for (let t = 0; t < tScroll; ++t) {
      for (let i of rowsInd) {
        rdKitCellRenderer.render(ctx, x, y, w, h, grid.cell(colName, i), null);
      }
    }
  });

  rowsInd = getIdxRandomSubarray(N, nScroll + tScroll);
  await testCase(`Vertical scrolling (${nScroll} random molecules, ${tScroll} times)`, () => {
    for (let t = 0; t < tScroll; ++t) {
      for (let i of rowsInd.slice(t, t + nScroll)) {
        rdKitCellRenderer.render(ctx, x, y, w, h, grid.cell(colName, i), null);
      }
    }
  });

  // Truncate the original dataset for searching
  rowInd = getIdxRandomSubset(N, nSearch);
  let newDf = DG.DataFrame.create();
  newDf.columns.addNew(colName, DG.TYPE.STRING);
  for (i of rowInd) {
    newDf.rows.addNew([col.get(i)]);
  }
  df = newDf;
  col = df.col(colName);

  const searchFor = [
    'c1ccccc1', // Benzene
    'O=C(C)Oc1ccccc1C(=O)O' // Aspirin
  ];

  await testCase(`Substructure search, building a library of ${col.length} molecules`, async () =>
    await grok.chem.substructureSearch(col));
  await testCase(`Substructure search, searching benzene in ${nSearch} molecules`, async () =>
    await grok.chem.substructureSearch(col, searchFor[0]));
  await testCase(`Substructure search, searching aspirin in ${nSearch} molecules`, async () =>
    await grok.chem.substructureSearch(col, searchFor[1]));

  await testCase(`Similarity scoring, building a library of ${col.length} molecules`, async () =>
    await grok.chem.similarityScoring(col));
  const queryIdx = getIdxRandomSubset(nSearch, nSample);
  await testCase(`Similarity scoring, search for ${queryIdx.length} samples in ${nSearch} molecules`, async () => {
    for (let i of queryIdx) {
      await grok.chem.similarityScoring(col, col.get(i));
    }
  });

})();