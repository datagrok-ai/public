import {before, category, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { _testDimensionalityReducer } from './dimensionality-reduce-utils';
import { createTableView, readDataframe } from './utils';
import { _package } from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';


category('chemSpace', async () => {

    let smallDf: DG.DataFrame;

    before(async () => {
      if (!chemCommonRdKit.moduleInitialized) {
        chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
        await chemCommonRdKit.initRdKitModuleLocal();
      }
        smallDf =  await readDataframe('sar-small.csv');
      });

  test('TSNE', async () => {
    await _testDimensionalityReducer(smallDf.col('smiles')!, 't-SNE');
  });
  test('UMAP', async () => {
    await _testDimensionalityReducer(smallDf.col('smiles')!, 'UMAP');
  });
  test('SPE', async () => {
    await _testDimensionalityReducer(smallDf.col('smiles')!, 'SPE');
  });

});
