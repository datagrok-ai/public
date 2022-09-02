import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { category, expect, test } from '@datagrok-libraries/utils/src/test';
import { _testSearchSubstructure, _testSearchSubstructureAllParameters } from './utils';
import { _testFindSimilar, _testGetSimilarities } from './menu-tests-similarity-diversity'
import { testCsv, testSubstructure } from './substructure-search-tests'

import {_importSdf} from '../open-chem/sdf-importer';

category('server features', () => {
  test('testDescriptors', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  });
});

category('chem exported', () => {
  
  test('findSimilar.api.sar-small', async () => {
    await _testFindSimilar(grok.chem.findSimilar);
  });

  test('getSimilarities.api.molecules', async () => {
    await _testGetSimilarities(grok.chem.getSimilarities);
  });

  test('substructureSearch', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    const trueIndices = [0, 2];
    const bitset: DG.BitSet = await grok.chem.searchSubstructure(df.col('smiles')!, testSubstructure);
    const bitsetString = bitset.toBinaryString();
    const bitsetArray = [...bitsetString];
    for (let k = 0; k < trueIndices.length; k++) {
      expect(bitsetArray[trueIndices[k]] === '1', true);
      bitsetArray[trueIndices[k]] = '0';
    }
  });

  test('mcs', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
      O=C1CN=C(c2ccccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
      O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    const mcs = await grok.chem.mcs(t.col('smiles')!);
    expect(mcs, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });
});