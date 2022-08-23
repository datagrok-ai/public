import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { category, expect, test } from '@datagrok-libraries/utils/src/test';
import { _testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText } from './utils';

import {_importSdf} from '../open-chem/sdf-importer';

category('chem', () => {

  test('testSubstructureSearch', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.searchSubstructure(t.col('smiles')!, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });

  test('testDescriptors', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  });

  test('testDiversitySearch', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.diversitySearch(t.col('smiles')!);
  });

  test('testMcs', async () => {
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

  test('chemblIdToSmiles', async () => {
    const query = 'Chembl:ChemblIdToSmiles';
    expect(await grok.functions.call(`${query}`, {'id': 'CHEMBL1185'}), 'CN(C)CCc1c[nH]c2ccc(C[C@H]3COC(=O)N3)cc12');
    expect(await grok.functions.call(`${query}`, {'id': 'CHEMBL45'}), 'COc1ccc2[nH]cc(CCNC(C)=O)c2c1');
    expect(await grok.functions.call(`${query}`, {'id': 'CHEMBL7784'}), 'CC(C)(CCCOc1ccc(OCCCC(C)(C)C(=O)O)c(-c2ccccc2Cl)c1)C(=O)O');
    expect(await grok.functions.call(`${query}`, {'id': 'CHEMBL6781'}), 'CCOC(=O)c1cc2ccccn2c(=O)n1');
    expect(await grok.functions.call(`${query}`, {'id': 'CHEMBL9812'}), 'COCCOC1(C2=NCCN2)COc2ccccc2O1');
  });
  
  test('molregnoToSmiles', async () => {
    const query = 'Chembl:molregnoToSmiles';
    expect(await grok.functions.call(`${query}`, {'molregno': 123}), 'O=c1oc2c(O)c(O)ccc2c2cc(F)ccc12');
    expect(await grok.functions.call(`${query}`, {'molregno': 241}), 'O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O');
    expect(await grok.functions.call(`${query}`, {'molregno': 1189}), 'CC(=O)Nc1nc(O)c2cc(S(=O)(=O)c3ccc4ccccc4c3)ccc2n1');
    expect(await grok.functions.call(`${query}`, {'molregno': 6190}), 'Cc1nn(CCCN2CCN(c3cccc(Cl)c3)CC2)c(=O)c2noc(C)c12');
    expect(await grok.functions.call(`${query}`, {'molregno': 5872}), 'N[C@@H](Cc1ccccc1)C(O)[C@@H](N)Cc1ccccc1');
  })
});
