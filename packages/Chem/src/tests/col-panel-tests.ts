import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';


const csvForInchi = `smiles
COc1ccc(CN(Cc2ccccc2)Cc3ccc(Br)cc3)cc1O
COc1ccc(CN(CCc2ccc(Br)cc2)Cc3ccccc3)cc1O
CCCCCCCC(=O)NCCC1CC(CC)(CC)C(=O)O1
CC1=C(C(C(=C(C)N1)C(=O)OCc2ccccc2)c3csc(n3)c4ccc(Cl)cc4)C(=O)
CCCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC)
CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OC(C)C)
CCOC(=O)C1=C(C)NC(=C(C1c2csc(n2)c3ccc(Cl)cc3)C(=O)OCC(C)C)`;

category('column panel', () => {
  test('inchi', async () => {
    if (DG.Test.isInBenchmark) {
      const df = await grok.data.files.openTable('System:AppData/Chem/tests/smiles_200K.zip');
      await grok.functions.call('Chem:addInchisTopMenu', {table: df, molecules: df.col('smiles')!});
    } else {
      await testInchiPanel('addInchisTopMenu', 'inchi', 'inchi',
        'InChI=1S/C22H22BrNO2/c1-26-22-12-9-19(13-21(22)25)16-24(14-17-5-3-2-4-6-17)15-18-7-10-20(23)11-8-18/h2-13,25H,14-16H2,1H3');
    }
  }, {benchmark: true});

  test('add inchi keys', async () => {
    if (DG.Test.isInBenchmark) {
      const df = await grok.data.files.openTable('System:AppData/Chem/tests/smiles_200K.zip');
      await grok.functions.call('Chem:addInchisKeysTopMenu', {table: df, molecules: df.col('smiles')!});
    } else
      await testInchiPanel('addInchisKeysTopMenu', 'inchiKeys', 'inchi_key', 'BRJLNESNMQYEGX-UHFFFAOYSA-N');
  }, {benchmark: true});
});

async function testInchiPanel(funcName: string, tableName: string, newColName: string, expected: string) {
  const t = DG.DataFrame.fromCsv(csvForInchi);
  t.name = tableName;
  const v = grok.shell.addTableView(t);
  await awaitCheck(() => v.dataFrame.name === tableName, undefined, 2000);
  await grok.functions.call(`Chem:${funcName}`, {table: t, molecules: t.columns.byName('smiles')});
  await delay(1000);
  expect(t.columns.names().includes(newColName), true, `${newColName} column has not been added`);
  expect(t.col(newColName)!.get(0) === expected, true);
  v.close();
}
