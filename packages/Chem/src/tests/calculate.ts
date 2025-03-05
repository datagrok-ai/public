import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, before, after, expect, test, delay, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {CONTAINER_TIMEOUT, ensureContainerRunning} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import { getMapIdentifiers } from '../widgets/identifiers';
import { addChemPropertiesColumns, addChemRisksColumns, addInchisKeysTopMenu, addInchisTopMenu, structuralAlertsTopMenu } from '../package';

category('calculate', () => {
  let smiles: DG.DataFrame;

  before(async () => {
        if (!chemCommonRdKit.moduleInitialized) {
          chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
          await chemCommonRdKit.initRdKitModuleLocal();
        }
    grok.shell.closeAll();
  });


  test('map identifiers', async () => {
    await ensureContainerRunning('name = "chem-chem"');
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);

    const checkIdentifier = async (name: string, value: string) => {
      await getMapIdentifiers(smiles, smiles.col('smiles')!, 'smiles', name);
      expect(smiles.columns.names().includes(name), true, `${name} column hasn't been added`);
      expect(smiles.col(name)?.get(2), value, `incorrect ${name} value`);
    }

    await checkIdentifier('inchi', 'InChI=1S/C6H4O3S2/c7-5-3-4(6(8)9-5)11-2-1-10-3/h1-2H2');
    await checkIdentifier('actor', '10489-75-5');

    const chemblPackInstalled = DG.Func.find({package: 'ChemblApi', name: 'getCompoundsIds'}).length;
    if (chemblPackInstalled) {
      await checkIdentifier('chembl', 'CHEMBL2262349');
      await checkIdentifier('pubchem', '82669'); 
      await checkIdentifier('mcule', 'MCULE-6494517749');     
    }
  }, {timeout: 60000 + CONTAINER_TIMEOUT});


  test('to inchi', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    addInchisTopMenu(smiles, smiles.col('smiles')!);
    expect(smiles.columns.names().includes('inchi'), true, `inchi column hasn't been added`);
    expect(smiles.col('inchi')?.get(5), 'InChI=1S/C13H9Cl3N4/c14-5-8-1-3-9(4-2-8)6-20-7-17-12-10(20)11(15)18-13(16)19-12/h1-4,7H,5-6H2', `incorrect inchi value`);
  }, {stressTest: true});

  test('to inchi keys', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    addInchisKeysTopMenu(smiles, smiles.col('smiles')!);
    expect(smiles.columns.names().includes('inchi_key'), true, `inchi_key column hasn't been added`);
    expect(smiles.col('inchi_key')?.get(5), 'BSKYPUUFGNPHLM-UHFFFAOYSA-N', `incorrect inchi_key value`);
  }, {stressTest: true});


  test('toxicity risks', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await addChemRisksColumns(smiles, smiles.col('smiles')!, true, true, true, true);
    expect(['Mutagenicity', 'Tumorigenicity', 'Irritating effects', 'Reproductive effects'].every((it) => smiles.columns.names().includes(it)),
      true, 'Not all toxicity columns have been added');
    expect(smiles.col('Mutagenicity')?.get(2), 'None', `incorrect Mutagenicity value`);
    expect(smiles.col('Tumorigenicity')?.get(2), 'None', `incorrect Tumorigenicity value`);
    expect(smiles.col('Irritating effects')?.get(2), 'High', `incorrect Irritating effects value`);
    expect(smiles.col('Reproductive effects')?.get(2), 'None', `incorrect Reproductive effects value`);
  }, {stressTest: true});


  test('properties', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await addChemPropertiesColumns(smiles, smiles.col('smiles')!, true, true, true);
    expect(['MW', 'HBA', 'HBD'].every((it) => smiles.columns.names().includes(it)),
      true, 'Not all properties columns have been added');
    expect(smiles.col('MW')?.get(3), 342.107177734375, `incorrect MW value`);
    expect(smiles.col('HBA')?.get(3), 5, `incorrect HBA value`);
    expect(smiles.col('HBD')?.get(3), 1, `incorrect HBD value`);
  });


  test('structural alerts', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await structuralAlertsTopMenu(smiles, smiles.col('smiles')!, true, true, true, false, false, false, false, false);
    expect(['PAINS (smiles)', 'BMS (smiles)', 'SureChEMBL (smiles)'].every((it) => smiles.columns.names().includes(it)),
      true, 'Not all structural alerts columns have been added');
    expect(smiles.col('PAINS (smiles)')?.get(2), false, `incorrect PAINS (smiles) value`);
    expect(smiles.col('BMS (smiles)')?.get(2), true, `incorrect BMS (smiles) value`);
    expect(smiles.col('SureChEMBL (smiles)')?.get(2), true, `incorrect SureChEMBL (smiles) value`);
  });


  after(async () => {
    grok.shell.closeAll();
  });
});


async function testGroup(groupName: string, funcName: string, colName: string, dlgName: string) {
  console.log(`testGroup, dlgName: ${dlgName}`);
  const smiles = grok.data.demo.molecules(20);
  const v = grok.shell.addTableView(smiles);
  await grok.data.detectSemanticTypes(smiles);
  await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
  grok.shell.topMenu.find('Chem').group(groupName).find(funcName).click();
  await getDlgAndClickOK(`cannot load ${funcName} dialog`, dlgName);
  await awaitCheck(() => v.dataFrame.columns.names().includes(colName), `${colName} column has not been added`, 10000);
  v.close();
  grok.shell.o = ui.div();
}


async function getDlgAndClickOK(error: string, header: string) {
  const dlg = () => {
    return Array.from(document.getElementsByClassName('d4-dialog'))
      .filter((dlg) => dlg.getElementsByClassName('d4-dialog-header')[0]
        .children[0].textContent?.toLocaleLowerCase() === header.toLocaleLowerCase());
  };
  await awaitCheck(() => {
    return dlg().length > 0;
  }, error, 5000);
  await delay(1000);
  Array.from(dlg()[0].getElementsByTagName('span')).find((el) => el.textContent === 'OK')?.click();
}
