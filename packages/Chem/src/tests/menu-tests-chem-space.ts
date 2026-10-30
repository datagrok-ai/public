import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package-test';
import {readDataframe} from './utils';
import {before, after, expect, category, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {chemSpace, runChemSpace} from '../analysis/chem-space';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getSimilaritiesMarix, getSimilaritiesMarixFromDistances} from '../utils/similarity-utils';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {MALFORMED_DATA_WARNING_CLASS} from '../constants';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import { getGPUDevice } from '@datagrok-libraries/math/src/webGPU/getGPUDevice';
import {PackageFunctions} from '../package';


category('top menu chem space', async () => {
  let smallDf: DG.DataFrame;
  let spgi100: DG.DataFrame;
  let approvedDrugs100: DG.DataFrame;
  let gd: GPUDevice | null;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    smallDf = await readDataframe('tests/sar-small_test.csv');
    spgi100 = await readDataframe('tests/spgi-100.csv');
    approvedDrugs100 = await readDataframe('tests/approved-drugs-100.csv');
    gd = await getGPUDevice();
  });

  test('chemSpaceOpens.smiles', async () => {
    const df = DG.Test.isInBenchmark ? gd ? await readDataframe('tests/smi10K.csv') : grok.data.demo.molecules(1000) : smallDf;
    await _testChemSpaceReturnsResult(df, 'smiles');
  }, {timeout: 1000000, benchmark: true, benchmarkTimeout: 1000000});

  test('chemSpaceOpens.molV2000', async () => {
    await _testChemSpaceReturnsResult(spgi100, 'Structure');
  });

  test('chemSpaceOpens.molV3000', async () => {
    await _testChemSpaceReturnsResult(approvedDrugs100, 'molecule');
  });

  test('chemSpace.emptyValues', async () => {
    const sarSmallEmptyRows = await readDataframe('tests/sar-small_empty_vals.csv');
    await _testChemSpaceReturnsResult(sarSmallEmptyRows, 'smiles');
  });

  test('chemSpace.malformedData', async () => {
    const testSmilesMalformed = await readDataframe('tests/Test_smiles_malformed.csv');
    DG.Balloon.closeAll();
    await _testChemSpaceReturnsResult(testSmilesMalformed, 'canonical_smiles');
    try {
      await waitUntil(() => {
        return document.querySelector(`.${MALFORMED_DATA_WARNING_CLASS}`)?.innerHTML ===
        '3 molecules with indexes 14,31,41 are possibly malformed and are not included in analysis';
      }, 'cannot find warning balloon', 5000);
    } finally {DG.Balloon.closeAll();}
  });

  test('chemSpace: the fingerprint options of the dialog reach the fingerprints', async () => {
    // the dialog's preprocessing settings (getFingerprints' inputs: type, radius, size, chirality, features)
    // a table of its own: a clone of smallDf would hold the caches the tests before filled
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', smallDf.col('smiles')!.toList())]);
    df.name = 'chem space options';
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    try {
      const sp = await runChemSpace(df, df.getCol('smiles'), DimReductionMethods.UMAP, BitArrayMetricsNames.Tanimoto,
        true, {preprocessingFuncArgs: {fingerprintType: 'Morgan', radius: 3, size: 1024}, fastRowCount: 100000});
      expect(sp != null, true);
      const caches = df.columns.names().filter((n) => n.startsWith('~smiles.Morgan'));
      expect(caches.length, 1, caches.join(', '));
      expect(caches[0].startsWith('~smiles.Morgan(radius=3, size=1024)'), true, caches[0]);
      expect(df.col(caches[0])!.get(0).length * 8, 1024, 'the fingerprints of 1024 bits');
    } finally {tv.close();}
  });

  test('chemSpace.allEmpty', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('id', Int32Array.from({length: 10}, (_, i) => i + 1)),
      DG.Column.fromStrings('structure', new Array(10).fill('')),
    ]);
    df.col('structure')!.semType = DG.SEMTYPE.MOLECULE;
    const tv = grok.shell.addTableView(df);
    try {
      const sp = await PackageFunctions.chemSpaceTopMenu(df, df.col('structure')!, DimReductionMethods.UMAP,
        BitArrayMetricsNames.Tanimoto, true, undefined, undefined, true);
      expect(sp?.type, DG.VIEWER.SCATTER_PLOT);
      for (const prefix of ['Embed_X_', 'Embed_Y_']) {
        const col = df.columns.toList().find((c) => c.name.startsWith(prefix));
        expect(col?.stats.missingValueCount, 0, `${prefix} column is missing or has empty values`);
      }
      expect(df.col('Cluster (DBSCAN)')?.length, 10);
      expect(df.rowCount, 10);
    } finally {tv.close();}
  });
});

async function _testChemSpaceReturnsResult(df: DG.DataFrame, col: string) {
  await grok.data.detectSemanticTypes(df);
  const tv = grok.shell.addTableView(df);
  await waitUntil(() => tv.name?.toLowerCase() === df.name?.toLowerCase(),
    'Chem space table view hasn\'t been created', 1000);
  try {
    const sp = await runChemSpace(df, df.getCol(col), DimReductionMethods.UMAP,
      BitArrayMetricsNames.Tanimoto, true, {fastRowCount: 100000});
    expect(sp != null, true);
  } finally {tv.close();}
}

