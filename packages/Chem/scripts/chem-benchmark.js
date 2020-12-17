//name: chemBenchmark
//language: javascript

const randomInt = (min, max) =>  min + Math.floor((max - min) * Math.random());
const range = (l, r) => new Array(r - l).fill().map((_,k) => k + l);

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
  let arr = range(0, N - 1);
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
  await time(caseCntr++ + '. ' + title, 1, f);
}

(async () => {

  const N_render = 1000;
  const N_scroll = 20;
  const t_scroll = 100;
  const N_search = 20000;
  const N_sample = 10;
  const x = 0, y = 0;
  const w = 200, h = 100;

  let df = await grok.data.getDemoTable('chem/zbb/99_p3_4.5-6.csv');
  const colName = 'smiles';
  let col = df.col(colName);
  let N = col.length;

  console.log('Chem Benchmark');

  let rdKitCellRenderer = new RDKitCellRenderer();
  let grid = DG.Viewer.grid(df);
  let canvas = grid.canvas;
  let ctx = canvas.getContext('2d');
  let rowsInd = null;

  rowsInd = getIdxRandomSubset(N, N_render);
  await testCase(`Rendering a ${rowsInd.length} random molecules`, () => {
    for (i of rowsInd) {
      rdKitCellRenderer.render(ctx, x, y, w, h, grid.cell(colName, i), null);
    }
  });

  rowsInd = getIdxRandomSubarray(N, N_scroll);
  await testCase(`Horizontal scrolling (${rowsInd.length} random molecules, ${t_scroll} times)`, () => {
    for (let t = 0; t < t_scroll; ++t) {
      for (let i of rowsInd) {
        rdKitCellRenderer.render(ctx, x, y, w, h, grid.cell(colName, i), null);
      }
    }
  });

  rowsInd = getIdxRandomSubarray(N, N_scroll + t_scroll);
  await testCase(`Vertical scrolling (${N_scroll} random molecules, ${t_scroll} times)`, () => {
    for (let t = 0; t < t_scroll; ++t) {
      for (let i of rowsInd.slice(t, t + N_scroll)) {
        rdKitCellRenderer.render(ctx, x, y, w, h, grid.cell(colName, i), null);
      }
    }
  });

  rowInd = getIdxRandomSubset(N, N_search);
  let newDf = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.STRING, colName, [])]);
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
    await grok.chem.substructureSearch(col, ''));
  await testCase(`Substructure search, searching benzene in ${N_search} molecules`, async () =>
    await grok.chem.substructureSearch(col, searchFor[0]));
  await testCase(`Substructure search, searching aspirin in ${N_search} molecules`, async () =>
    await grok.chem.substructureSearch(col, searchFor[1]));

  await testCase(`Similarity scoring, building a library of ${col.length} molecules`, async () =>
      await grok.chem.similarityScoring(col, ''));
  const queryIdx = getIdxRandomSubset(N, N_sample);
  await testCase(`Similarity scoring, search for ${queryIdx.length} samples in ${N_search} molecules`, async () => {
    for (let i of queryIdx) {
      await grok.chem.similarityScoring(col, col.get(i));
    }
  });

})();