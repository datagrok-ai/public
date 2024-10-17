import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package-test';
import {readDataframe} from './utils';
import {before, after, expect, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {chemSpace, runChemSpace} from '../analysis/chem-space';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getSimilaritiesMarix, getSimilaritiesMarixFromDistances} from '../utils/similarity-utils';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {MALFORMED_DATA_WARNING_CLASS} from '../constants';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import { getGPUDevice } from '@datagrok-libraries/math/src/webGPU/getGPUDevice';


category('top menu chem space', async () => {
  let smallDf: DG.DataFrame;
  let spgi100: DG.DataFrame;
  let approvedDrugs100: DG.DataFrame;
  let gd = await getGPUDevice();

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    smallDf = await readDataframe('tests/sar-small_test.csv');
    spgi100 = await readDataframe('tests/spgi-100.csv');
    approvedDrugs100 = await readDataframe('tests/approved-drugs-100.csv');
  });

  test('chemSpaceOpens.smiles', async () => {
    const df = DG.Test.isInBenchmark ? gd ? await grok.data.files
      .openTable('System:AppData/Chem/tests/smiles_100K.zip') : await readDataframe('tests/smiles_50K.csv') : smallDf;
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
      await awaitCheck(() => {
        return document.querySelector(`.${MALFORMED_DATA_WARNING_CLASS}`)?.innerHTML ===
        '2 molecules with indexes 31,41 are possibly malformed and are not included in analysis';
      },
      'cannot find warning balloon', 5000);
    } finally {DG.Balloon.closeAll();}
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
});

async function _testChemSpaceReturnsResult(df: DG.DataFrame, col: string) {
  await grok.data.detectSemanticTypes(df);
  const tv = grok.shell.addTableView(df);
  await awaitCheck(() => tv.name?.toLowerCase() === df.name?.toLowerCase(),
    'Chem space table view hasn\'t been created', 1000);
  try {
    const sp = await runChemSpace(df, df.getCol(col), DimReductionMethods.UMAP,
      BitArrayMetricsNames.Tanimoto, true, {fastRowCount: 100000});
    expect(sp != null, true);
  } finally {tv.close();}
}

