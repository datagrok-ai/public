import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { category, expect, test } from '@datagrok-libraries/utils/src/test';
import { _testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText } from './utils';
import { _testFindSimilar, _testGetSimilarities } from './menu-tests-similarity-diversity'

import {_importSdf} from '../open-chem/sdf-importer';

category('server features', () => {

  test('findSimilarServer.api.sar-small', async () => {
    await _testFindSimilar(grok.chem.findSimilarServer);
  });

  // test('testSubstructureSearch', async () => {
  //   const t = grok.data.demo.molecules();
  //   await grok.chem.searchSubstructure(t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  // });

  // test('testDescriptors', async () => {
  //   const t = grok.data.demo.molecules();
  //   await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  // });

  // test('testDiversitySearch', async () => {
  //   const t = grok.data.demo.molecules();
  //   await grok.chem.diversitySearch(t.col('smiles')!);
  // });

  test('mcs', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
O=C1CN=C(c2ccccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    await grok.chem.mcs(t.col('smiles')!);
  });
});

category('chem exported', () => {
  
  test('findSimilar.api.sar-small', async () => {
    await _testFindSimilar(grok.chem.findSimilar);
  });

  test('getSimilarities.api.molecules', async () => {
    await _testGetSimilarities(grok.chem.getSimilarities);
  });
});