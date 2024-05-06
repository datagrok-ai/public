import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import { awaitCheck, before, category, delay, expect, test } from '@datagrok-libraries/utils/src/test';
import { createTableView } from './utils';
import { mmpAnalysis } from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import { _package } from '../package-test';
import { MmpAnalysis } from '../analysis/molecular-matched-pairs/mmp-analysis';

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
        'Mean Difference Activity': {idxs: [0, 11, 30], values: [-2.5343997478485107, 3.5564699172973633, 3.458324432373047]},
        'Mean Difference Permeability': {idxs: [0, 11, 30], values: [1.8472713232040405, -4.707592010498047, -7.125835418701172]},
        'Mean Difference Toxicity': {idxs: [0, 11, 30], values: [1.4948477745056152, -0.617708146572113, -0.6957154273986816]},
    },
    'Transformations_Pairs': {
        'From': {idxs: [0, 30, 50], values: [
            pairsFromMolblock,
            'O=c1cc(-CC(C)C)oc2cc(O)c(O)c(O)c12',
            'O=c1cc(-CC)oc2cc(O)c(O)c(O)c12'
            ]},
        'To': {idxs: [0, 30, 50], values: [
            pairsToMolblock, 
            'O=c1cc(-CCC(=O)NC)oc2cc(O)c(O)c(O)c12', 
            'O=c1cc(-C(C)Br)oc2cc(O)c(O)c(O)c12'
        ]},
        'Difference Activity': {idxs: [0, 30, 50], values: [-3.014911651611328, 4.207612991333008, 8.298847198486328]},
        'Difference Permeability': {idxs: [0, 30, 50], values: [2.2153091430664062, -12.706818580627441, -8.750051498413086]},
        'Difference Toxicity': {idxs: [0, 30, 50], values: [1.5138638019561768, 0.28581464290618896, 0.2180633544921875]},
    }
}

category('mmpa', () => {

    before(async () => {
        grok.shell.closeAll();
        if (!chemCommonRdKit.moduleInitialized) {
            chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
            await chemCommonRdKit.initRdKitModuleLocal();
        }
    });

    test('mmpa opens', async () => {
        const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
        await mmpAnalysis(tv.dataFrame, tv.dataFrame.col('smiles')!,
            tv.dataFrame.clone().columns.remove('smiles'));
        //ensure MMPA opened
        await awaitCheck(() => document.getElementsByClassName('chem-mmpa-transformation-tab-header').length > 0,
            'MMPA hasn\'t been started', 3000);
        //ensure fragments and pairs grids have been created
        await awaitCheck(() => document.getElementsByClassName('d4-grid').length === 3,
            'Fragments and Pairs grids haven\'t been created', 3000);
        //ensure embeddings columns have been created for cliffs tab
        await awaitCheck(() => tv.dataFrame.columns.names().includes('~Embed_X_1')
            && tv.dataFrame.columns.names().includes('~Embed_Y_1'), 'Embeddings haven\'t been created', 3000);
        //ensure embeddings columns have been calculated
        await awaitCheck(() => tv.dataFrame.col('~Embed_X_1')!.stats.missingValueCount === 0
            && tv.dataFrame.col('~Embed_Y_1')!.stats.missingValueCount === 0 , 'Embeddings haven\'t been calculated', 5000);
    });

    test('transformations tab', async () => {
        const tv = await createTableView('demo_files/matched_molecular_pairs.csv');
        const mmp: MmpAnalysis = await mmpAnalysis(tv.dataFrame, tv.dataFrame.col('smiles')!,
            tv.dataFrame.clone().columns.remove('smiles'));
        await delay(5000);

        //check Fragments Grid
        const fragsDf = mmp.allPairsGrid.dataFrame;
        await awaitCheck(() => fragsDf.rowCount === 40 && fragsDf.columns.length === 7
            && fragsDf.filter.trueCount === 2 && fragsDf.filter.get(0) && fragsDf.filter.get(2),
            'Incorrect fragments grid', 3000);
        checkRandomValues(fragsDf, 'Transformations_Fragments');

        //check Pairs Grid
        const pairsDf = mmp.casesGrid.dataFrame;
        await awaitCheck(() => pairsDf.rowCount === 54 && pairsDf.columns.length === 13
            && pairsDf.filter.trueCount === 3 && pairsDf.filter.get(0) && pairsDf.filter.get(1) && pairsDf.filter.get(2),
            'Incorrect pairs grid', 3000);
        checkRandomValues(mmp.casesGrid.dataFrame, 'Transformations_Pairs');

        //changing fragment
        mmp.allPairsGrid.dataFrame.currentRowIdx = 2;
        await awaitCheck(() => pairsDf.filter.trueCount === 2 && pairsDf.filter.get(6) && pairsDf.filter.get(7),
            'Pairs haven\'t been changed after fragment change', 3000);

        //changing target molecule
        tv.dataFrame.currentRowIdx = 4;
        await awaitCheck(() => fragsDf.filter.trueCount === 3
            && fragsDf.filter.get(3) && fragsDf.filter.get(4) && fragsDf.filter.get(7)
            && pairsDf.filter.trueCount === 2 && pairsDf.filter.get(8) && pairsDf.filter.get(9),
            'Pairs haven\'t been changed after fragment change', 3000);
    });

});


function checkRandomValues(df: DG.DataFrame, dfName: string) {
    Object.keys(randomValsToCheck[dfName]).forEach((key: string) => {
        const idxs = randomValsToCheck[dfName][key].idxs;
        const vals = randomValsToCheck[dfName][key].values;
        idxs.forEach((it, idx) => expect(df.col(key)!.get(it), vals[idx], `incorrect data in ${key} column, row ${it}`));
    })
}