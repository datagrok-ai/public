import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {awaitCheck, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {createTableView} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import {MatchedMolecularPairsViewer} from '../analysis/molecular-matched-pairs/mmp-viewer/mmp-viewer';

const pairsFromMolblock = `
     RDKit          2D

 13 14  0  0  0  0  0  0  0  0999 V2000
    0.6347   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1346   -0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6347    0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3653    0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8654    1.2990    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8654   -0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3653   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8655   -2.1650    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3654   -3.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6346   -3.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1346   -2.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1346   -2.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6346   -3.0311    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
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
    0.6346   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1345   -0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6347    0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3653    0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8654    1.2990    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8653   -0.4330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3654   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8653   -2.1650    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3654   -3.0312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6346   -3.0312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1345   -2.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1345   -2.1650    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6346   -3.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6346   -3.0310    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
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
    'From': {idxs: [1, 29, 37], values: ['CC[*:1]', 'CNC(=O)C[*:1]', 'CC(Br)[*:1]']},
    'To': {idxs: [5, 9, 37], values: ['O[*:1]', 'Br[*:1]', 'CC[*:1]']},
    'Pairs': {idxs: [0, 10, 39], values: [3, 2, 1]},
    'Mean Difference Activity':
      {idxs: [0, 11, 30], values: [-2.53, 3.56, 3.46]},
    'Mean Difference Permeability':
      {idxs: [0, 11, 30], values: [1.85, -4.71, -7.13]},
    'Mean Difference Toxicity':
      {idxs: [0, 11, 30], values: [1.49, -0.62, -0.70]},
  },
  'Transformations_Pairs': {
    'From': {idxs: [0, 30, 50], values: [
      pairsFromMolblock,
      'O=c1cc(-CC(C)C)oc2cc(O)c(O)c(O)c12',
      'O=c1cc(-CC)oc2cc(O)c(O)c(O)c12',
    ]},
    'To': {idxs: [0, 30, 50], values: [
      pairsToMolblock,
      'O=c1cc(-CCC(=O)NC)oc2cc(O)c(O)c(O)c12',
      'O=c1cc(-C(C)Br)oc2cc(O)c(O)c(O)c12',
    ]},
    'Difference Activity':
      {idxs: [0, 30, 50], values: [-3.01, 4.21, 8.30]},
    'Difference Permeability':
      {idxs: [0, 30, 50], values: [2.22, -12.71, -8.75]},
    'Difference Toxicity': {idxs: [0, 30, 50], values: [1.51, 0.29, 0.22]},
  },
  'Generation': {
    'Structure': {idxs: [0, 70, 113], values: [
      'c12ccc(Cl)cc2nccc1CC',
      'O=C1Oc2ccccc2C(O)C1Cc3ccc(O)cc3O',
      'c12ccc(Cl)cc2nccc1CO',
    ]},
    'Initial value': {idxs: [0, 70, 113], values: [48.67, 5.89, -2.46]},
    'Activity': {idxs: [0, 38, 76], values: ['Activity', 'Permeability', 'Toxicity']},
    'Core': {idxs: [0, 70, 113],
      values: ['Clc1ccc2c(C[*:1])ccnc2c1', 'O=C1Oc2ccccc2C([*:1])C1Cc1ccc(O)cc1O', 'Clc1ccc2c(C[*:1])ccnc2c1']},
    'From': {idxs: [0, 70, 107], values: ['C[*:1]', 'O[*:1]', 'O=C(O)[*:1]']},
    'To': {idxs: [0, 70, 113], values: ['O=C(O)[*:1]', 'CC(C)[*:1]', 'CC[*:1]']},
    'Prediction': {idxs: [0, 70, 113], values: [58.98, 13.53, 0.33]},
    'Generation': {idxs: [0, 70, 113],
      values: ['O=C(O)Cc1ccnc2cc(Cl)ccc12', 'CC(C)C1c2ccccc2OC(=O)C1Cc1ccc(O)cc1O', 'CCCc1ccnc2cc(Cl)ccc12']},
  },
};


const randomValsToCheckGPU: {[key: string]: {[key: string]: {idxs: number[], values: any[]}}} = {
  'Transformations_Fragments': {
    'From': {idxs: [1, 9, 29, 37], values: ['C[*:1]', 'Br[*:1]', 'CNC(=O)C[*:1]', 'CC(Br)[*:1]']},
    'To': {idxs: [1, 5, 9, 29, 37], values: ['CC[*:1]', 'O[*:1]', 'Cl[*:1]', 'C[*:1]', 'CC[*:1]']},
    'Pairs': {idxs: [0, 10, 39], values: [3, 2, 1]},
    'Mean Difference Activity':
      {idxs: [0, 11, 30], values: [2.53, 3.56, 3.46]},
    'Mean Difference Permeability':
      {idxs: [0, 11, 30], values: [1.85, -4.71, -7.13]},
    'Mean Difference Toxicity':
      {idxs: [0, 11, 30], values: [1.49, -0.62, -0.70]},
  },
  'Transformations_Pairs': {
    'From': {idxs: [3, 30, 50], values: [
      pairsFromMolblock,
      'O=c1cc(-CC(C)C)oc2cc(O)c(O)c(O)c12',
      'O=c1cc(-CC)oc2cc(O)c(O)c(O)c12',
    ]},
    'To': {idxs: [3, 30, 50], values: [
      pairsToMolblock,
      'O=c1cc(-CCC(=O)NC)oc2cc(O)c(O)c(O)c12',
      'O=c1cc(-C(C)Br)oc2cc(O)c(O)c(O)c12',
    ]},
    'Difference Activity': {idxs: [0, 30, 50], values: [3.01, 4.21, 8.30]},
    'Difference Permeability':
      {idxs: [0, 30, 50], values: [-2.22, -12.71, -8.75]},
    'Difference Toxicity':
      {idxs: [0, 30, 50], values: [-1.51, 0.29, 0.22]},
  },
  'Generation': {
    'Structure': {idxs: [0, 70, 113], values: [
      'c12ccc(Cl)cc2nccc1CC',
      'O=C1Oc2ccccc2C(O)C1Cc3ccc(O)cc3O',
      'c12ccc(Cl)cc2nccc1CO',
    ]},
    'Initial value': {idxs: [0, 70, 113], values: [48.67, 5.89, -2.46]},
    'Activity': {idxs: [0, 38, 76], values: ['Activity', 'Permeability', 'Toxicity']},
    'Core': {idxs: [0, 70, 113],
      values: ['Clc1ccc2c(C[*:1])ccnc2c1', 'O=C1Oc2ccccc2C([*:1])C1Cc1ccc(O)cc1O', 'Clc1ccc2c(C[*:1])ccnc2c1']},
    'From': {idxs: [0, 70, 107], values: ['C[*:1]', 'O[*:1]', 'O=C(O)[*:1]']},
    'To': {idxs: [0, 70, 113], values: ['O=C(O)[*:1]', 'CC(C)[*:1]', 'CC[*:1]']},
    'Prediction': {idxs: [0, 70, 113], values: [58.98, 13.53, 0.33]},
    'Generation': {idxs: [0, 70, 113],
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
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
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
    //ensure MMPA opened
    await awaitCheck(() => document.getElementsByClassName('chem-mmpa-transformation-tab-header').length > 0,
      'MMPA hasn\'t been started', 3000);
    //ensure fragments and pairs grids have been created
    await awaitCheck(() => document.getElementsByClassName('d4-grid').length === 3,
      'Fragments and Pairs grids haven\'t been created', 3000);
    //ensure embeddings columns have been created for cliffs tab
    await awaitCheck(() => tv.dataFrame.columns.names().includes('~Embed_X_1') &&
      tv.dataFrame.columns.names().includes('~Embed_Y_1'), 'Embeddings haven\'t been created', 3000);
    //ensure embeddings columns have been calculated
    await awaitCheck(() => tv.dataFrame.col('~Embed_X_1')!.stats.missingValueCount === 0 &&
      tv.dataFrame.col('~Embed_Y_1')!.stats.missingValueCount === 0, 'Embeddings haven\'t been calculated', 10000);
    expect(mmp.mmpa!.rules.rules.length, 40, `Incorrect rules`);
    expect(mmp.mmpa!.rules!.smilesFrags.length, 14, `Incorrect smilesFrags`);
  });

  test('transformationsTab', async () => {
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
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

    //check Fragments Grid
    await awaitCheck(() => mmp.pairedGrids!.fpGrid.dataFrame != null, 'All pairs grid hasn\'t been created', 3000);
    const fragsDf = mmp.pairedGrids!.fpGrid.dataFrame;
    await awaitCheck(() => fragsDf.rowCount === 40 && fragsDf.columns.length === 7 &&
      fragsDf.filter.trueCount === 2 && fragsDf.filter.get(mmp.calculatedOnGPU ? 1 : 0) && fragsDf.filter.get(2),
    'Incorrect fragments grid', 3000);
    checkRandomValues(fragsDf, 'Transformations_Fragments', mmp.calculatedOnGPU!);

    //check Pairs Grid
    const pairsDf = mmp.pairedGrids!.mmpGrid.dataFrame;
    await awaitCheck(() => pairsDf.rowCount === 54 && pairsDf.columns.length === 13 &&
    pairsDf.filter.trueCount === 3 && pairsDf!.filter.get(mmp.calculatedOnGPU ? 3 : 0) &&
    pairsDf.filter.get(mmp.calculatedOnGPU ? 4 : 1) && pairsDf.filter.get(mmp.calculatedOnGPU ? 5 : 2),
    'Incorrect pairs grid', 3000);
    checkRandomValues(mmp.pairedGrids!.mmpGrid.dataFrame, 'Transformations_Pairs', mmp.calculatedOnGPU!);

    //changing fragment
    mmp.pairedGrids!.fpGrid.dataFrame.currentRowIdx = 2;
    await awaitCheck(() => pairsDf.filter.trueCount === 2 && pairsDf.filter.get(6) && pairsDf.filter.get(7),
      'Pairs haven\'t been changed after fragment change', 3000);

    //changing target molecule
    tv.dataFrame.currentRowIdx = 4;
    await awaitCheck(() => fragsDf.filter.trueCount === 3 &&
        fragsDf.filter.get(3) && fragsDf.filter.get(4) && fragsDf.filter.get(7) &&
        pairsDf.filter.trueCount === 2 && pairsDf.filter.get(8) && pairsDf.filter.get(9),
    'Pairs haven\'t been changed after fragment change', 3000);
  });

  test('cliffsTab', async () => {
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
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
    //check created lines
    await awaitCheck(() => mmp.lines!.from.length === 81 && mmp.lines!.to.length === 81 &&
    mmp.linesIdxs!.length === 81, 'Incorrect lines number');
    await awaitCheck(() => mmp.linesMask!.allTrue, 'Incorrect initial lines mask');
    checkRandomArrayVals(mmp.lines!.from, [0, 10, 30, 50, 70], [30, 6, 37, 23, 9], 'mmp.lines.from');
    checkRandomArrayVals(mmp.lines!.to, [0, 10, 30, 50, 70], [0, 28, 0, 27, 23], 'mmp.lines.to');
    checkRandomArrayVals(mmp.linesIdxs!,
      [0, 10, 30, 50, 80], [mmp.calculatedOnGPU ? 0 : 3, 22, 8, 47, 52], 'mmp.linesIdxs');
    checkRandomArrayVals(mmp.lines!.colors,
      [0, 30, 80], ['31,119,180', '255,187,120', '44,160,44'], 'mmp.lines.colors');
    checkRandomArrayVals(mmp.linesActivityCorrespondance, [0, 27, 55], [0, 1, 2], 'mmp.linesActivityCorrespondance');

    //changing sliders inputs values
    mmp.mmpFilters!.activitySliderInputs![0].value = 11.87;
    mmp.mmpFilters!.activitySliderInputs![1].value = 14.15;
    mmp.mmpFilters!.activitySliderInputs![2].value = 1.627;
    await awaitCheck(() => DG.BitSet.fromBytes(mmp.linesMask!.buffer.buffer, 81).trueCount === 7,
      'Incorrect lines mask after slider input changed', 3000);

    //switch of one of activities
    mmp.mmpFilters!.activityActiveInputs[2].value = false;
    await awaitCheck(() => DG.BitSet.fromBytes(mmp.linesMask!.buffer.buffer, 81).trueCount === 2,
      'Incorrect lines mask after checkboxes values changed', 3000);
  });

  test('generationTab', async () => {
    const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
    let mmp: MatchedMolecularPairsViewer = (grok.shell.v as DG.TableView)
      .addViewer('Matched Molecular Pairs Analysis') as MatchedMolecularPairsViewer;
    mmp.setOptions({
      molecules: 'smiles',
      activities: tv.dataFrame.clone().columns.remove('smiles').names(),
      fragmentCutoff: 0.4});
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
    await awaitCheck(() => mmp.generationsGrid?.dataFrame != null, 'Generation grid hasn\'t been created', 3000);
    const genDf = mmp.generationsGrid!.dataFrame;
    await awaitCheck(() => genDf.rowCount === 114, 'Incorrect lines number');
    checkRandomValues(genDf, 'Generation', mmp.calculatedOnGPU!);
    //check that 'Existing' column has been calculated
    await awaitCheck(() =>
      genDf.columns.names().includes('Existing'), '\'Existing\' column hasn\'t been created', 10000);

    const isExpected = genDf.col('Existing')!.toList().filter((it) => it).length == 22 ||
      genDf.col('Existing')!.toList().filter((it) => it).length == 23;
    expect(isExpected, true, 'Incorrect data in \'Existing\' column');
  });
});


function checkRandomValues(df: DG.DataFrame, dfName: string, isGPU: boolean) {
  const randVals = isGPU ? randomValsToCheckGPU : randomValsToCheck;
  Object.keys(randVals[dfName]).forEach((key: string) => {
    const idxs = randVals[dfName][key].idxs;
    const vals = randVals[dfName][key].values;
    idxs.forEach((it, idx) => {
      const val = typeof df.col(key)!.get(it) === 'number' ?
        Math.round(df.col(key)!.get(it)*100)/100 :
        df.col(key)!.get(it);
      expect(val, vals[idx], `incorrect data in ${key} column, row ${it}`);
    });
  });
}

function checkRandomArrayVals(array: any, idxs: number[], vals: (number | string)[], name: string) {
  idxs.forEach((it: number, idx: number) => expect(array[it], vals[idx], `Incorrect value in ${name}, idx: ${it}`));
}
