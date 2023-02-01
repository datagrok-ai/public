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

  test('add inchi panel', async () => {
    await testInchiPanel('addInchisTopMenu', 'inchi');
  });

  test('inchi keys', async () => {
    await testInchiPanel('addInchisKeysTopMenu', 'inchiKeys');
  });
});

async function testInchiPanel(funcName: string, tableName: string) {
  const t = DG.DataFrame.fromCsv(csvForInchi);
  t.name = tableName;
  const v = grok.shell.addTableView(t);
  await awaitCheck(() => v.name === tableName, undefined, 1000);
  await grok.functions.call(`Chem:${funcName}`, {table: t, molecules: t.columns.byName('smiles')});
  v.close();
}
