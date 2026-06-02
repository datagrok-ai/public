import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {awaitCheck, before, category, expect, test} from './_timed-test';
import {createTableView} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import {MatchedMolecularPairsViewer} from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer';
import { MMP_NAMES, SHOW_FRAGS_MODE } from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-constants';

const _tsLog = (msg: string): void => console.log(`[${new Date().toISOString()}] ${msg}`);

const pairsFromMolblock = `
     RDKit          2D

 13 14  0  0  0  0  0  0  0  0999 V2000
    0.4546    1.4171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9545    2.2831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4546    3.1491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5454    3.1491    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0455    4.0151    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.0455    2.2831    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5454    1.4171    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0455    0.5511    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5454   -0.3149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4546   -0.3149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9545    0.5511    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9545    0.5511    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4546   -0.3149    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  4  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  1  0
 11  1  1  0
  7  1  1  0
M  END
`;

const pairsToMolblock = `
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
    0.4546    1.4172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9545    2.2832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4546    3.1492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5454    3.1492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0455    4.0152    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.0453    2.2832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5454    1.4172    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0454    0.5512    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5455   -0.3150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4545   -0.3150    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9544    0.5512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9544    0.5512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4545   -0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4545   -0.3148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  4  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 11  1  1  0
  7  1  1  0
M  END
`;

const randomValsToCheck: {[key: string]: {[key: string]: {idxs: number[], values: any[]}}} = {
  'Transformations_Fragments': {
    'From': {idxs: [1, 29, 37], values: ["CC(Br)[*:1]", "C[*:1]", "O[*:1]"]},
    'To': {idxs: [5, 9, 37], values: ["CC(C)C[*:1]", "CC(C)[*:1]", "O[*:1]"]},
    'Total count': {idxs: [0, 28, 39], values: [1, 2, 3]},
    '\u0394 Activity':
      {idxs: [0, 11, 30], values: [-12.23, -4.64, 6.85]},
    '\u0394 Permeability':
      {idxs: [0, 11, 30], values: [-14.59, -7.56, 8.35]},
    '\u0394 Toxicity':
      {idxs: [0, 11, 30], values: [-2.79, -0.65, 0.86]},
  },
  'Transformations_Pairs': {
    'From': {idxs: [0, 30, 50], values: [
      "CC(Br)C1C(=O)Oc2ccccc2C1O",
      "CCCc1cc(=O)c2c(O)c(O)c(O)cc2o1",
      "O=C1Oc2ccccc2C(O)C1CO",
    ]},
    'To': {idxs: [0, 30, 50], values: [
      "CC(Br)C1C(=O)Oc2ccccc2C1O",
      "CCCc1cc(=O)c2c(O)c(O)c(O)cc2o1",
      "O=C1Oc2ccccc2C(O)C1CO",
    ]},
    '\u0394 Activity':
      {idxs: [0, 30, 50], values: [-12.23, 1.78, 10.31]},
    '\u0394 Permeability':
      {idxs: [0, 30, 50], values: [-14.59, 1.23, 12.71]},
    '\u0394 Toxicity': {idxs: [0, 30, 50], values: [-3.19, 0.29, 1.68]},
  },
  'Generation': {
    'Structure': {idxs: [0, 70, 109], values: [
      'CCc1ccnc2cc(Cl)ccc12',
      'O=C1Oc2ccccc2C(O)C1Cc1ccc(O)cc1O',
      'OCc1ccnc2cc(Cl)ccc12',
    ]},
    [MMP_NAMES.OBSERVED]: {idxs: [0, 70, 109], values: [48.67, 5.89, -2.46]},
    [MMP_NAMES.PROPERTY_TYPE]: {idxs: [0, 38, 76], values: ['Activity', 'Permeability', 'Toxicity']},
    [MMP_NAMES.CORE]: {idxs: [0, 70, 109],
      values: ['Clc1ccc2c(C[*:1])ccnc2c1', 'O=C1Oc2ccccc2C([*:1])C1Cc1ccc(O)cc1O', 'Clc1ccc2c(C[*:1])ccnc2c1']},
    [MMP_NAMES.FROM]: {idxs: [0, 70, 107], values: ['C[*:1]', 'O[*:1]', 'O[*:1]']},
    [MMP_NAMES.TO]: {idxs: [0, 70, 109], values: ['O=C(O)[*:1]', 'CC(C)[*:1]', 'CC[*:1]']},
    [MMP_NAMES.PREDICTED]: {idxs: [0, 70, 109], values: [58.98, 13.53, 0.33]},
    [MMP_NAMES.NEW_MOLECULE]: {idxs: [0, 70, 109],
      values: ['O=C(O)Cc1ccnc2cc(Cl)ccc12', 'CC(C)C1c2ccccc2OC(=O)C1Cc1ccc(O)cc1O', 'CCCc1ccnc2cc(Cl)ccc12']},
  },
};


category('mmpa', () => {
  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('mmpaOpens', async () => {
    _tsLog('[MMPA] mmpaOpens: loading table view');
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    _tsLog(`[MMPA] mmpaOpens: table loaded, rows=${tv.dataFrame.rowCount}, cols=${tv.dataFrame.columns.length}`);
    _tsLog('[MMPA] mmpaOpens: adding MMPA viewer');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    _tsLog('[MMPA] mmpaOpens: viewer added, calling setOptions');
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
    _tsLog('[MMPA] mmpaOpens: setOptions returned, awaiting MMPA in tv.viewers');
    //wait for MMPA viewer to be created
    await awaitCheck(() => {
      let mmpCreated = false;
      for (const v of tv.viewers) {
        if (v.type === 'Matched Molecular Pairs Analysis') {
          mmp = v as MatchedMolecularPairsViewer;
          mmpCreated = true;
        }
      }
      return mmpCreated;
    }, 'MMPA hasn\'t been initialized', 3000);
    _tsLog('[MMPA] mmpaOpens: MMPA viewer present in tv.viewers, awaiting transformation tab');
    //ensure MMPA opened
    await awaitCheck(() => document.getElementsByClassName('chem-mmpa-transformation-tab-header').length > 0,
      'MMPA hasn\'t been started', 3000);
    _tsLog('[MMPA] mmpaOpens: transformation tab present, awaiting 3 d4-grid elements');
    //ensure fragments and pairs grids have been created
    await awaitCheck(() => document.getElementsByClassName('d4-grid').length === 3,
      'Fragments and Pairs grids haven\'t been created', 3000);
    _tsLog(`[MMPA] mmpaOpens: grids present, mmpa=${!!mmp.mmpa}, ` +
      `rules=${mmp.mmpa?.rules?.rules?.length}, smilesFrags=${mmp.mmpa?.rules?.smilesFrags?.length}`);
    expect(mmp.mmpa!.rules.rules.length, 40, `Incorrect rules`);
    expect(mmp.mmpa!.rules!.smilesFrags.length, 14, `Incorrect smilesFrags`);
    _tsLog('[MMPA] mmpaOpens: done');
  });

  test('transformationsTab', async () => {
    // dispose the previous test's view: a leftover live MMPA/cliffs scatterplot
    // keeps re-rendering and starves the event loop, freezing the next test's render
    grok.shell.closeAll();
    _tsLog('[MMPA] transformationsTab: loading table view');
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    _tsLog('[MMPA] transformationsTab: adding MMPA viewer');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    _tsLog('[MMPA] transformationsTab: calling setOptions');
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
    _tsLog('[MMPA] transformationsTab: awaiting MMPA in tv.viewers');
    //wait for MMPA viewer to be created
    await awaitCheck(() => {
      let mmpCreated = false;
      for (const v of tv.viewers) {
        if (v.type === 'Matched Molecular Pairs Analysis') {
          mmp = v as MatchedMolecularPairsViewer;
          mmpCreated = true;
        }
      }
      return mmpCreated;
    }, 'MMPA hasn\'t been initialized', 3000);

    _tsLog('[MMPA] transformationsTab: awaiting fpGrid.dataFrame');
    //check Fragments Grid
    await awaitCheck(() => mmp.pairedGrids?.fpGrid?.dataFrame != null, 'All pairs grid hasn\'t been created', 3000);
    const fragsDf = mmp.pairedGrids!.fpGrid.dataFrame;
    _tsLog(`[MMPA] transformationsTab: fragsDf ready, awaiting rowCount=40,cols=10 ` +
      `(got ${fragsDf.rowCount}/${fragsDf.columns.length})`);
    await awaitCheck(() => fragsDf.rowCount === 40 && fragsDf.columns.length === 10, 'Incorrect fragments grid', 3000);
    _tsLog('[MMPA] transformationsTab: checking Transformations_Fragments values');
    checkRandomValues(fragsDf, 'Transformations_Fragments', true);

    _tsLog('[MMPA] transformationsTab: awaiting mmpGridTrans.dataFrame');
    //check Pairs Grid
    await awaitCheck(() => mmp.pairedGrids?.mmpGridTrans?.dataFrame != null, 'mmpGrid hasn\'t been created', 3000);
    const pairsDf = mmp.pairedGrids!.mmpGridTrans.dataFrame;
    _tsLog(`[MMPA] transformationsTab: pairsDf ready, awaiting rowCount=54,cols=19 ` +
      `(got ${pairsDf.rowCount}/${pairsDf.columns.length})`);
    await awaitCheck(() => pairsDf.rowCount === 54 && pairsDf.columns.length === 19, 'Incorrect pairs grid', 3000);
    _tsLog('[MMPA] transformationsTab: checking Transformations_Pairs values');
    checkRandomValues(mmp.pairedGrids!.mmpGridTrans.dataFrame, 'Transformations_Pairs', true);

    _tsLog('[MMPA] transformationsTab: changing fragGrid current row to 2');
    //changing fragment
    mmp.followCurrentRowInFragGrid!.value = true;
    mmp.pairedGrids!.fpGrid.dataFrame.currentRowIdx = 2;
    await awaitCheck(() => pairsDf.filter.trueCount === 2 && pairsDf.filter.get(6) && pairsDf.filter.get(7),
      'Pairs haven\'t been changed after fragment change', 3000);
    _tsLog('[MMPA] transformationsTab: pairs filter updated after fragment change');

    mmp.showFragmentsChoice!.value = SHOW_FRAGS_MODE.Current;
    _tsLog('[MMPA] transformationsTab: changing target molecule row to 4');
    //changing target molecule
    tv.dataFrame.currentRowIdx = 4;
    await awaitCheck(() => fragsDf.filter.trueCount === 3 &&
        fragsDf.filter.get(3) && fragsDf.filter.get(4) && fragsDf.filter.get(7),
    'Fragments haven\'t been changed after target change', 3000);
    _tsLog('[MMPA] transformationsTab: done');
  });

  test('cliffsTab', async () => {
    grok.shell.closeAll();
    _tsLog('[MMPA] cliffsTab: loading table view');
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    _tsLog('[MMPA] cliffsTab: adding MMPA viewer');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    _tsLog('[MMPA] cliffsTab: calling setOptions');
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
    _tsLog('[MMPA] cliffsTab: awaiting MMPA in tv.viewers');
    //wait for MMPA viewer to be created
    await awaitCheck(() => {
      let mmpCreated = false;
      for (const v of tv.viewers) {
        if (v.type === 'Matched Molecular Pairs Analysis') {
          mmp = v as MatchedMolecularPairsViewer;
          mmpCreated = true;
        }
      }
      return mmpCreated;
    }, 'MMPA hasn\'t been initialized', 3000);
    _tsLog('[MMPA] cliffsTab: awaiting Cliffs tab DOM element (10s timeout)');
    //open cliffs tab
    await awaitCheck(() => mmp.root.querySelector('[name=\'Cliffs\']') != null,
      'Cliffs tab hasn\'t been created', 10000);
    _tsLog('[MMPA] cliffsTab: clicking Cliffs tab');
    const cliffsTabHeader = mmp.root.querySelector('[name=\'Cliffs\']') as HTMLElement;
    cliffsTabHeader.click();
    _tsLog('[MMPA] cliffsTab: awaiting ~Embed_X_1 and ~Embed_Y_1 columns');
    //ensure embeddings columns have been created for cliffs tab
    await awaitCheck(() => tv.dataFrame.columns.names().includes('~Embed_X_1') &&
     tv.dataFrame.columns.names().includes('~Embed_Y_1'), 'Embeddings haven\'t been created', 3000);
    _tsLog('[MMPA] cliffsTab: awaiting embeddings calculated (10s timeout)');
    //ensure embeddings columns have been calculated
    await awaitCheck(() => tv.dataFrame.col('~Embed_X_1')!.stats.missingValueCount === 0 &&
     tv.dataFrame.col('~Embed_Y_1')!.stats.missingValueCount === 0, 'Embeddings haven\'t been calculated', 10000);
    _tsLog('[MMPA] cliffsTab: awaiting lines (expecting 81)');
    //check created lines
    await awaitCheck(() => mmp.lines?.from.length === 81 && mmp.lines?.to.length === 81 &&
    mmp.linesIdxs!.length === 81, 'Incorrect lines number', 3000);
    _tsLog('[MMPA] cliffsTab: awaiting initial lines mask');
    await awaitCheck(() => mmp.linesMask?.allTrue == false, 'Incorrect initial lines mask');
    _tsLog('[MMPA] cliffsTab: checking line array values');
    checkRandomArrayVals(mmp.lines?.from, [0, 10, 30, 50, 70], [30, 6, 37, 23, 9], 'mmp.lines.from');
    checkRandomArrayVals(mmp.lines?.to, [0, 10, 30, 50, 70], [0, 28, 0, 27, 23], 'mmp.lines.to');
    checkRandomArrayVals(mmp.linesIdxs,
      [0, 10, 30, 50, 80], [mmp.calculatedOnGPU ? 0 : 3, 22, 8, 47, 52], 'mmp.linesIdxs');
    checkRandomArrayVals(mmp.lines?.colors,
      [0, 30, 80], ['31,119,180', '255,187,120', '44,160,44'], 'mmp.lines.colors');
    checkRandomArrayVals(mmp.linesActivityCorrespondance, [0, 27, 55], [0, 1, 2], 'mmp.linesActivityCorrespondance');

    expect(mmp.mmpFilters?.activitySliderInputs?.length, 3, 'mmp cliffs filters haven\'t been created');
    _tsLog('[MMPA] cliffsTab: changing activity slider values');
    //changing sliders inputs values
    mmp.mmpFilters!.activitySliderInputs![0].value = 11.87;
    mmp.mmpFilters!.activitySliderInputs![1].value = 14.15;
    mmp.mmpFilters!.activitySliderInputs![2].value = 1.627;
    await awaitCheck(() => !!mmp.linesMask && DG.BitSet.fromBytes(mmp.linesMask.buffer.buffer, 81).trueCount === 7,
      'Incorrect lines mask after slider input changed', 3000);
    _tsLog('[MMPA] cliffsTab: lines mask updated after slider change, toggling activity[2] off');

    //switch of one of activities
    mmp.mmpFilters!.activityActiveInputs[2].value = false;
    await awaitCheck(() => !!mmp.linesMask && DG.BitSet.fromBytes(mmp.linesMask!.buffer.buffer, 81).trueCount === 2,
      'Incorrect lines mask after checkboxes values changed', 3000);
    _tsLog('[MMPA] cliffsTab: done');
  });

  test('generationTab', async () => {
    grok.shell.closeAll();
    _tsLog('[MMPA] generationTab: loading table view');
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    _tsLog('[MMPA] generationTab: adding MMPA viewer');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    _tsLog('[MMPA] generationTab: calling setOptions');
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
    _tsLog('[MMPA] generationTab: awaiting MMPA in tv.viewers');
    //wait for MMPA viewer to be created
    await awaitCheck(() => {
      let mmpCreated = false;
      for (const v of tv.viewers) {
        if (v.type === 'Matched Molecular Pairs Analysis') {
          mmp = v as MatchedMolecularPairsViewer;
          mmpCreated = true;
        }
      }
      return mmpCreated;
    }, 'MMPA hasn\'t been initialized', 3000);
    _tsLog('[MMPA] generationTab: awaiting Generation tab DOM element (10s timeout)');
    await awaitCheck(() => mmp.root.querySelector('[name=\'Generation\']') != null,
      'Generation tab hasn\'t been created', 10000);
    _tsLog('[MMPA] generationTab: clicking Generation tab');
    const genTabHeader = mmp.root.querySelector('[name=\'Generation\']') as HTMLElement;
    genTabHeader.click();
    _tsLog('[MMPA] generationTab: awaiting generationsGrid.dataFrame (10s timeout)');
    await awaitCheck(() => mmp.generationsGrid?.dataFrame != null, 'Generation grid hasn\'t been created', 10000);
    const genDf = mmp.generationsGrid!.dataFrame;
    _tsLog(`[MMPA] generationTab: genDf present, awaiting rowCount=110 (got ${genDf.rowCount})`);
    await awaitCheck(() => genDf.rowCount === 110, 'Incorrect row count');
    _tsLog('[MMPA] generationTab: checking Generation values');
    checkRandomValues(genDf, 'Generation');
    _tsLog('[MMPA] generationTab: awaiting Prediction column (10s timeout)');
    //check that 'Prediction' column has been calculated
    await awaitCheck(() =>
      genDf.columns.names().includes('Prediction'), '\'Prediction\' column hasn\'t been created', 10000);

    const isExpected = () => {
      const predicted = genDf.col('Prediction')!.toList().filter((it) => it).length;
      // we get different results on cpu and gpu for predictions
      return predicted === 88 || (mmp.calculatedOnGPU && predicted === 89);
    };
    const predictedCount = genDf.col('Prediction')!.toList().filter((it) => it).length;
    _tsLog(`[MMPA] generationTab: Prediction filled count=${predictedCount}, ` +
      `calculatedOnGPU=${mmp.calculatedOnGPU}`);
    expect(isExpected(), true, 'Incorrect data in \'Prediction\' column');
    _tsLog('[MMPA] generationTab: done');
  });
});


function checkRandomValues(df: DG.DataFrame, dfName: string, sort?: boolean) {
  Object.keys(randomValsToCheck[dfName]).forEach((key: string) => {
    const idxs = randomValsToCheck[dfName][key].idxs;
    const vals = randomValsToCheck[dfName][key].values;
    const col = df.col(key)?.toList();
    if (sort)
      typeof col![0] === 'number' ? col?.sort((a, b) => a - b) : col?.sort();
    idxs.forEach((it, idx) => {
      const val = typeof col![it] === 'number' ?
        Math.round(col![it]*100)/100 : col![it];
      expect(val, vals[idx], `incorrect data in ${key} column, row ${it}`);
    });
  });
}

function checkRandomArrayVals(array: any, idxs: number[], vals: (number | string)[], name: string) {
  expect(array != undefined, true, 'array of values is not defined');
  idxs.forEach((it: number, idx: number) => expect(array[it], vals[idx], `Incorrect value in ${name}, idx: ${it}`));
}
