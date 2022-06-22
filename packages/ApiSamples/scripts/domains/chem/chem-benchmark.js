//name: chemBenchmark
//language: javascript

// https://stackoverflow.com/a/47593316
function mulberry32(a) {
    return function() {
      var t = a += 0x6D2B79F5;
      t = Math.imul(t ^ t >>> 15, t | 1);
      t ^= t + Math.imul(t ^ t >>> 7, t | 61);
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    }
}
var rand = mulberry32(42); // seed = 42

const randomInt = (min, max) => min + Math.floor((max - min) * rand());
const range = (l, r) => new Array(r - l).fill().map((_, k) => k + l);

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
  let msg = await DG.timeAsync(++caseCntr + '. ' + title, f);
  if (typeof msg === 'string')
  	console.log("   " + msg);
}

(async () => {

  const nSearch = 25000;
  const nRender = 1000;
  const nScroll = 20;
  const tScroll = 100;
  const nSample = 10;

  let df = await grok.data.getDemoTable('chem/chembl/chembl-100k.csv');
  // For using with your HOME files in Datagrok:
  // let df = (await grok.functions.eval('OpenServerFile("UserName:Home/Chembl_100K.csv")'))[0];
  // An alternative data source for testing (~40K+ molecules): chem/zbb/99_p3_4.5-6.csv
  if (nSearch < df.rowCount)
    df.rows.removeAt(nSearch, df.rowCount - nSearch, false);
  const colName = 'smiles';
  let col = df.col(colName);

  console.clear();
  console.log('Chem Benchmark');

  let rdKitCellRenderer = await grok.functions.call('Chem:rdKitCellRenderer');
  let grid = DG.Viewer.grid(df);
  let canvas = grid.canvas;
  let ctx = canvas.getContext('2d');
  let rowsInd = null;

  const render = (i) => rdKitCellRenderer.render(
    ctx, 0, 0, 200, 100, grid.cell(colName, i), null);

  rowsInd = getIdxRandomSubarray(nSearch, nRender);
  await testCase(`Rendering a ${rowsInd.length} random molecules`, () => {
    for (i of rowsInd) {
      render(i);
    }
  });

  rowsInd = getIdxRandomSubarray(nSearch, nScroll);
  await testCase(`Horizontal scrolling (${rowsInd.length} random molecules, ${tScroll} times)`, () => {
    for (let t = 0; t < tScroll; ++t) {
      for (let i of rowsInd) {
        render(i);
      }
    }
  });

  rowsInd = getIdxRandomSubarray(nSearch, nScroll + tScroll);
  await testCase(`Vertical scrolling (${nScroll} random molecules, ${tScroll} times)`, () => {
    for (let t = 0; t < tScroll; ++t) {
      for (let i of rowsInd.slice(t, t + nScroll)) {
        render(i);
      }
    }
  });

  const searchFor = [
    'c1ccccc1', // Benzene https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL277500/ 
    'CC(=O)Oc1ccccc1C(=O)O', // Aspirin https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25/
    'CCCCCCCCCCCCCCCCCC(=O)O', // Stearic acid https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL46403/
    'NCCCC[C@H](N)C(=O)O' // Lusine https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL8085/
  ];
  
  const searchForNames = ['Benzene', 'Aspirin', 'Stearic acid', 'Lusine'];

  await testCase(`Substructure search, building a library of ${col.length} molecules`, async () =>
    await grok.chem.searchSubstructure(col));
  
  
  for (let i = 0; i < searchFor.length; ++i)
    await testCase(`Substructure search, searching ${searchForNames[i]} in ${nSearch} molecules`, async () => {
      let s = await grok.chem.searchSubstructure(col, searchFor[i]); return `Found ${s.trueCount} molecules for ${searchFor[i]}`; });

  await testCase(`Similarity scoring, building a library of ${col.length} molecules`, async () =>
    await grok.chem.getSimilarities(col));
  const queryIdx = getIdxRandomSubset(nSearch, nSample);
  await testCase(`Similarity scoring, search for ${queryIdx.length} samples in ${nSearch} molecules`, async () => {
    for (let i of queryIdx) {
      await grok.chem.getSimilarities(col, col.get(i));
    }
  });
  
  await testCase(`Substructure search (server), searching benzene in ${nSearch} molecules`, async () =>
    await grok.chem.searchSubstructureServer(col, searchFor[0]), false);
  await testCase(`Substructure search (server), searching aspirin in ${nSearch} molecules`, async () =>
    await grok.chem.searchSubstructureServer(col, searchFor[1]), false);
  
  /*
  await testCase(`Similarity scoring (server), search for ${queryIdx.length} samples in ${nSearch} molecules`, async () => {
    await grok.chem.findSimilarServer(col, searchFor[0]);
  });
  */


})();