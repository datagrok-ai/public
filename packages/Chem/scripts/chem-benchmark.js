//name: chemBenchmark
//language: javascript

const randomInt = (min, max) => min + Math.floor((max - min) * Math.random());
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
  await DG.timeAsync(++caseCntr + '. ' + title, f);
}

(async () => {

  const nSearch = 25000;
  const nRender = 1000;
  const nScroll = 20;
  const tScroll = 100;
  const nSample = 10;

  let df = await grok.data.getDemoTable('chem/chembl/chembl_100k.csv');
  // For using with your HOME files in Datagrok:
  // let df = (await grok.functions.eval('OpenServerFile("UserName:Home/Chembl_100K.csv")'))[0];
  // An alternative data source for testing (~40K+ molecules): chem/zbb/99_p3_4.5-6.csv
  if (nSearch < df.rowCount)
    df.rows.removeAt(nSearch, df.rowCount - nSearch, false);
  const colName = 'smiles';
  let col = df.col(colName);

  console.log('Chem Benchmark');

  // TODO: await _ChemPackage.init(); dummy priming for now
  await grok.chem.searchSubstructure(
    DG.Column.fromList(DG.TYPE.STRING, '', ['c1ccccc1']));
  let rdKitCellRenderer = new RDKitCellRenderer();
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
    'CC(=O)Oc1ccccc1C(=O)O' // Aspirin https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL25/
  ];

  await testCase(`Substructure search, building a library of ${col.length} molecules`, async () =>
    await grok.chem.searchSubstructure(col));
  await testCase(`Substructure search, searching benzene in ${nSearch} molecules`, async () =>
    await grok.chem.searchSubstructure(col, searchFor[0]));
  await testCase(`Substructure search, searching aspirin in ${nSearch} molecules`, async () =>
    await grok.chem.searchSubstructure(col, searchFor[1]));

  await testCase(`Similarity scoring, building a library of ${col.length} molecules`, async () =>
    await grok.chem.getSimilarities(col));
  const queryIdx = getIdxRandomSubset(nSearch, nSample);
  await testCase(`Similarity scoring, search for ${queryIdx.length} samples in ${nSearch} molecules`, async () => {
    for (let i of queryIdx) {
      await grok.chem.getSimilarities(col, col.get(i));
    }
  });

})();