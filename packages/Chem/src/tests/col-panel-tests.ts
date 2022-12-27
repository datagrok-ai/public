import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {awaitCheck, category, delay, expect, expectFloat, test} from '@datagrok-libraries/utils/src/test';
import {_testSearchSubstructure,
  _testSearchSubstructureAllParameters,
  _testSearchSubstructureSARSmall,
  loadFileAsText} from './utils';

import {_importSdf} from '../open-chem/sdf-importer';

const csvForInchi = `smiles
COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O
COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O
CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1
CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)
CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)
CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)
CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)`

category('column panel', () => {
  test('mcs', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
      O=C1CN=C(c2ccccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
      O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    const v = grok.shell.addTableView(t);
    await grok.functions.call('Chem:addMcsPanel', {col: t.columns.byName('smiles')});
    v.close();
  });

  test('add inchi panel', async () => {
    await testInchiPanel('addInchisPanel', 'inchi');
  });

  test('inchi keys', async () => {
    await testInchiPanel('addInchisKeysPanel', 'inchiKeys');
  });
});

async function testInchiPanel(funcName: string, tableName: string) {
  const t = DG.DataFrame.fromCsv(csvForInchi);
  t.name = tableName;
  const v = grok.shell.addTableView(t);
  await awaitCheck(() => v.name === tableName, undefined, 1000);
  await grok.functions.call(`Chem:${funcName}`, {col: t.columns.byName('smiles')});
  v.close();
}
