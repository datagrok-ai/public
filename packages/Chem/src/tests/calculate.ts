import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, before, after, expect, expectArray, expectFloat, test, delay, awaitCheck}
  from '@datagrok-libraries/test/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/test/src/test-container-utils';
import {CONTAINER_TIMEOUT, readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';
import {getMapIdentifiers} from '../widgets/identifiers';
import {PackageFunctions} from '../package';
import {ChemicalEngine, cruxReady, setChemicalEngine} from '../engine/chem-engine';
import {getPropertiesAsColumns} from '../widgets/properties';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';

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
    await ensureContainerRunning('name = "chem-chem"', CONTAINER_TIMEOUT);
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);

    const checkIdentifier = async (name: string, value: string) => {
      await getMapIdentifiers(smiles, smiles.col('smiles')!, 'smiles', name);
      expect(smiles.columns.names().includes(name), true, `${name} column hasn't been added`);
      expect(smiles.col(name)?.get(2), value, `incorrect ${name} value`);
    };

    await checkIdentifier('inchi', 'InChI=1S/C6H4O3S2/c7-5-3-4(6(8)9-5)11-2-1-10-3/h1-2H2');

    const chemblPackInstalled = DG.Func.find({package: 'ChemblApi', name: 'getCompoundsIds'}).length;
    if (chemblPackInstalled) {
      await checkIdentifier('actor', '10489-75-5');
      await checkIdentifier('chembl', 'CHEMBL2262349');
      await checkIdentifier('pubchem', '82669');
      await checkIdentifier('mcule', 'MCULE-6494517749');
    }
  }, {timeout: 60000 + CONTAINER_TIMEOUT});


  test('to inchi', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await PackageFunctions.addInchisTopMenu(smiles, smiles.col('smiles')!);
    expect(smiles.columns.names().includes('inchi'), true, `inchi column hasn't been added`);
    expect(smiles.col('inchi')?.get(5), 'InChI=1S/C13H9Cl3N4/c14-5-8-1-3-9(4-2-8)6-20-7-17-12-10(20)11(15)18-13(16)19-12/h1-4,7H,5-6H2', `incorrect inchi value`);
  }, {stressTest: true});

  test('to inchi keys', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await PackageFunctions.addInchisKeysTopMenu(smiles, smiles.col('smiles')!);
    expect(smiles.columns.names().includes('inchi_key'), true, `inchi_key column hasn't been added`);
    expect(smiles.col('inchi_key')?.get(5), 'BSKYPUUFGNPHLM-UHFFFAOYSA-N', `incorrect inchi_key value`);
  }, {stressTest: true});


  test('toxicity risks', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await PackageFunctions.addChemRisksColumns(smiles, smiles.col('smiles')!, true, true, true, true);
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
    setChemicalEngine(ChemicalEngine.RDKit); // OpenChemLib's values (Crux mode's: the next test)
    try {
      await PackageFunctions.addChemPropertiesColumns(smiles, smiles.col('smiles')!, true, true, true);
    } finally {
      setChemicalEngine(null);
    }
    expect(['MW', 'HBA', 'HBD'].every((it) => smiles.columns.names().includes(it)),
      true, 'Not all properties columns have been added');
    expect(smiles.col('MW')?.get(3), 342.107177734375, `incorrect MW value`);
    expect(smiles.col('HBA')?.get(3), 5, `incorrect HBA value`);
    expect(smiles.col('HBD')?.get(3), 1, `incorrect HBD value`);
  });

  test('properties: Crux mode, RDKit\'s definitions', async () => {
    // the pinned RDKit's values (crux-core's tools/chem-wasm/chem_props_oracle.py): MolWt, not OpenChemLib's
    // monoisotopic mass; ESOL, Crippen, RDKit's strict rotatable bonds; empty cells where no molecule is read
    const col = DG.Column.fromStrings('smiles', ['CN1CCC(CC1)Oc2ccc(cc2)C(F)(F)F',
      'CC1CCN(CC(O)CN2C(=O)C3=C(SCCS3)C2=O)CC1', 'NCc1cnn(O)c1Br', '', 'C1CC(']);
    col.semType = DG.SEMTYPE.MOLECULE;
    const names = ['MW', 'HBA', 'HBD', 'LogP', 'LogS', 'PSA', 'Rotatable bonds', 'Stereo centers', 'Molecule charge',
      'Molecular formula'];
    const expected: {[name: string]: any[]} = {
      'MW': [259.271, 342.486, 192.016], 'HBA': [2, 5, 4], 'HBD': [0, 1, 2], 'LogP': [3.1784, 1.1396, 0.3416],
      'LogS': [-3.9460940582866364, -2.8353577331083555, -1.4954699830328941], 'PSA': [12.47, 111.45, 64.07],
      'Rotatable bonds': [2, 4, 1], 'Stereo centers': [0, 1, 0], 'Molecule charge': [0, 0, 0],
      'Molecular formula': ['C13H16F3NO', 'C15H22N2O3S2', 'C4H6BrN3O']};
    setChemicalEngine(ChemicalEngine.Crux);
    try {
      expect(await cruxReady(), true, 'the Crux module loads');
      const cols = await getPropertiesAsColumns(col, names);
      expect(cols.map((c) => `${c.name}:${c.type}`).join(', '), 'MW:double, HBA:int, HBD:int, LogP:double, ' +
        'LogS:double, PSA:double, Rotatable bonds:int, Stereo centers:int, Molecule charge:int, ' +
        'Molecular formula:string');
      for (const c of cols) {
        expected[c.name].forEach((value, i) => typeof value === 'number' && c.type === DG.TYPE.FLOAT ?
          expectFloat(c.get(i), value, 1e-4, `${c.name}, row ${i}`) : expect(c.get(i), value, `${c.name}, row ${i}`));
        expect(c.isNone(3) && c.isNone(4), true, `${c.name}: empty for an empty cell and a molecule not read`);
      }
    } finally {
      setChemicalEngine(null);
    }
  });


  test('structural alerts', async () => {
    smiles = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(smiles);
    await PackageFunctions.structuralAlertsTopMenu(smiles, smiles.col('smiles')!, true, true, true, false, false, false, false, false);
    expect(['PAINS (smiles)', 'BMS (smiles)', 'SureChEMBL (smiles)'].every((it) => smiles.columns.names().includes(it)),
      true, 'Not all structural alerts columns have been added');
    expect(smiles.col('PAINS (smiles)')?.get(2), false, `incorrect PAINS (smiles) value`);
    expect(smiles.col('BMS (smiles)')?.get(2), true, `incorrect BMS (smiles) value`);
    expect(smiles.col('SureChEMBL (smiles)')?.get(2), true, `incorrect SureChEMBL (smiles) value`);
  });
});


category('clustering', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('bitbirch', async () => {
    const df = await readDataframe('smiles.csv');
    await grok.data.detectSemanticTypes(df);
    await _bitbirch(df, df.col('canonical_smiles')!);
    const clusters = df.col('Cluster (BitBIRCH)')!;
    expect(clusters.stats.missingValueCount, 0);
    expect(clusters.stats.uniqueCount, 366);
    expect(df.rowCount, 1000);
  });

  test('similarity matrix', async () => {
    const df = await readDataframe('smiles.csv');
    df.rows.removeAt(50, df.rowCount - 50);
    await grok.data.detectSemanticTypes(df);
    const molecules = df.col('canonical_smiles')!;
    const matrix = await PackageFunctions.similarityMatrixTopMenu(df, molecules, df.col('molregno')!);
    try {
      expect(matrix.name, 'canonical_smiles similarity matrix');
      expect(matrix.rowCount, 50);
      expectArray(matrix.columns.names(), ['symbol', ...df.col('molregno')!.toList().map((v) => `${v}`)]);
      const cols = matrix.columns.toList().slice(1);
      let offDiagonal = 0;
      for (let i = 0; i < 50; i++) {
        expect(cols[i].get(i), 1, `diagonal cell ${i}`);
        for (let j = 0; j < 50; j++) {
          expectFloat(cols[j].get(i), cols[i].get(j), 1e-6, `cells ${i},${j} and ${j},${i}`);
          if (cols[j].get(i) < 1)
            offDiagonal++;
        }
      }
      expect(offDiagonal > 0, true, 'every pair of molecules is fully similar');
      const scores = (await PackageFunctions.getSimilarities(molecules, molecules.get(0))).columns.byIndex(0);
      for (let i = 1; i < 50; i++)
        expectFloat(cols[0].get(i), scores.get(i), 1e-4, `similarity of row ${i} to row 0`);
    } finally {
      grok.shell.tableView(matrix.name)?.close();
    }
  });

  test('cluster MCS', async () => {
    const df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
    const molecules = df.col('Structure')!;
    await _bitbirch(df, molecules);
    const clusters = df.col('Cluster (BitBIRCH)')!;
    const mcs = await PackageFunctions.performClusterMCS(molecules, clusters);
    expect(mcs.semType, DG.SEMTYPE.MOLECULE);
    expect(mcs.stats.missingValueCount, 0);
    expect(mcs.stats.uniqueCount, 81);
    const rdkit = PackageFunctions.getRdKitModule();
    for (let i = 0; i < df.rowCount; i++) {
      const mol = rdkit.get_mol(molecules.get(i));
      const query = getQueryMolSafe(mcs.get(i)!, '', rdkit)!;
      try {
        expect(mol.get_substruct_match(query) !== '{}', true,
          `row ${i} does not contain the MCS of its cluster, ${mcs.get(i)}`);
      } finally {
        mol.delete();
        query.delete();
      }
    }
  });
});

// through the registered function: the worker loads its wasm from a web root only the package bundle knows
async function _bitbirch(df: DG.DataFrame, molecules: DG.Column): Promise<void> {
  await DG.Func.find({package: 'Chem', name: 'bitbirchClusteringTopMenu'})[0]
    .prepare({table: df, molecules: molecules, threshold: 0.55, fingerprintType: 'Morgan'})
    .call(undefined, undefined, {processed: false});
}

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
