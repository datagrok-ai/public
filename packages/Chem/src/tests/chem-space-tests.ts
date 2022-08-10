import { before, category, test } from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import { _testChemSpaceReturnsResult, _testDimensionalityReducer } from './dimensionality-reduce-utils';
import { readDataframe } from './utils';
import { _package } from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';


category('chemSpace', async () => {

  let smallDf: DG.DataFrame;
  let sarSmallEmptyRows: DG.DataFrame;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    smallDf = await readDataframe('sar-small.csv');
    sarSmallEmptyRows = await readDataframe('tests/sar-small_empty_vals.csv');
  });


  test('chemSpaceOpens', async () => {
    await _testChemSpaceReturnsResult(smallDf, 'UMAP');
  });
  test('chemSpaceOpensWithEmptyRows', async () => {
    await _testChemSpaceReturnsResult(sarSmallEmptyRows, 'UMAP');
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
