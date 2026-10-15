import {category, expectArray, test, before, awaitCheck, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {PackageFunctions} from '../package';
import {_convertMolNotation} from '../utils/convert-notation-utils';

// import {_package} from '../package-test';
import {loadFileAsText, readDataframe} from './utils';

category('converters', async () => {
  let molfileV2K: string;
  let molfileV3K: string;
  let molecules: { [key: string]: string[]};
  let rdkitModule: any;
  before(async () => { // wait until RdKit module gets initialized
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    molfileV2K = await loadFileAsText('tests/molfileV2000_test_convert.mol');
    molfileV3K = await loadFileAsText('tests/molfileV3000_test_convert.mol');
    molecules = {
      smiles: ['CN1C(=O)CN=C(C2CCCCC2)c2ccccc21'],
      smarts: [
        '[#6]-[#7]1-[#6](=[#8])-[#6]-[#7]=[#6](-[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]-2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1',
      ],
      molblock: [molfileV2K],
      v3Kmolblock: [molfileV3K],
    };
    rdkitModule = await grok.functions.call('Chem:getRdKitModule');
  });


  function _testConvert(srcNotation: DG.chem.Notation, tgtNotation: DG.chem.Notation) {
    const result = [];
    for (const mol of molecules[srcNotation])
      result.push(_convertMolNotation(mol, srcNotation, tgtNotation, PackageFunctions.getRdKitModule()));
    expectArray(result.map((it) => it.replaceAll('\r', '')),
      molecules[tgtNotation].map((it) => it.replaceAll('\r', '')));
  }

  test('SMILES to Molfile V2000', async () => {
    if (DG.Test.isInBenchmark) {
      const df = await readDataframe('smiles.csv');
      const rdkitModule = PackageFunctions.getRdKitModule();
      for (let i = 0; i < df.rowCount; i++)
        _convertMolNotation(df.get('canonical_smiles', i), DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, rdkitModule);
    } else
      _testConvert(DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
  }, {benchmark: true});

  test('SMILES to SMARTS', async () => {
    _testConvert(DG.chem.Notation.Smiles, DG.chem.Notation.Smarts);
  });
  test('Molfile V2000 to SMILES', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
  });
  test('Molfile V2000 to SMARTS', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
  });
  test('Molfile V2000 to V3000', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock);
  });
  test('Molfile V3000 to SMARTS', async () => {
    _testConvert(DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smarts);
  });
  test('Molfile V3000 to SMILES', async () => {
    _testConvert(DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles);
  });
  test('Names to SMILES', async () => {
    const chemblPackInstalled = DG.Func.find({package: 'ChemblApi', name: 'getCompoundsIds'}).length;
    if (chemblPackInstalled) {
      const df = await readDataframe('tests/names_to_smiles.csv');
      await PackageFunctions.namesToSmiles(df, df.col('Name')!);
      await awaitCheck(() => df.columns.names().includes(`canonical_smiles`),
        `Column with names has not been converted to smiles`, 5000);
      expect(df.get('canonical_smiles', 4) === '');
      let mol1; let mol2; let smiles1; let smiles2;
      try {
        mol1 = rdkitModule.get_mol(df.get('canonical_smiles', 3));
        smiles1 = mol1.get_smiles();
        mol2 = rdkitModule.get_mol('O=C(CCCN1CCC(n2c(O)nc3ccccc32)CC1)c1ccc(F)cc1');
        smiles2 = mol1.get_smiles();
      } finally {
        mol1?.delete();
        mol2?.delete();
      }
      expect(smiles1 === smiles2);
    }
  });
  test('Convert notation round trip', async () => {
    const df = await readDataframe('smiles.csv');
    await grok.data.detectSemanticTypes(df);
    const smiles = df.col('canonical_smiles')!;
    await PackageFunctions.convertNotation(df, smiles, DG.chem.Notation.MolBlock);
    const molblocks = df.col('canonical_smiles_molblock')!;
    expect(molblocks.semType, DG.SEMTYPE.MOLECULE);
    expect(molblocks.toList().every((m) => m.includes('V2000') && m.includes('M  END')), true);
    await PackageFunctions.convertNotation(df, molblocks, DG.chem.Notation.Smiles);
    await PackageFunctions.recalculateCoords(df, smiles, 'CoordGen');
    const recalculated = df.col('canonical_smiles_recalcCoords')!;
    expect(recalculated.toList().every((m, i) => m.includes('M  END') && m !== molblocks.get(i)), true,
      'CoordGen kept the coordinates of a converted molblock');

    const canonical = (col: DG.Column, flat: boolean) => col.toList().map((s) => {
      const mol = rdkitModule.get_mol(s);
      const bare = flat ? rdkitModule.get_mol(mol.get_smiles().replace(/[@\/\\]/g, '')) : null;
      try {
        return (bare ?? mol).get_smiles();
      } finally {
        mol.delete();
        bare?.delete();
      }
    });
    // GROK-20956: RDKit 2024.09's molblock round trip inverts the vinyl stereocentre of rows 477, 486 and 487
    const molblockFlips = [477, 486, 487];
    // RDKit's CoordGen flips the C=N geometry of row 832
    const coordGenFlips = [832];
    const converted: [DG.Column, number[]][] = [[molblocks, molblockFlips],
      [df.col('canonical_smiles_molblock_smiles')!, molblockFlips], [recalculated, coordGenFlips]];
    for (const flat of [true, false]) {
      const original = canonical(smiles, flat);
      for (const [col, flips] of converted) {
        const values = canonical(col, flat);
        const differ = values.map((_, i) => i).filter((i) => values[i] !== original[i] && (flat || !flips.includes(i)));
        expect(differ.length, 0, `${col.name} holds other molecules${flat ? ' ignoring stereochemistry' : ''} ` +
          `in rows ${differ.slice(0, 10).join(', ')}`);
      }
    }
  }, {timeout: 60000});

  test('Convert notations for column', async () => {
    const df = DG.Test.isInBenchmark ? await readDataframe('tests/smi10K.csv') :
      await readDataframe('tests/Test_smiles_with_empty_and_malformed.csv');

    const testColumnAdded = async (colName: string, convertFrom: DG.chem.Notation, convertTo: DG.chem.Notation) => {
      await PackageFunctions.convertNotation(df, df.col(colName)!, convertTo, false);
      await awaitCheck(() => df.columns.names().includes(`${colName}_${convertTo}`),
        `Column with ${convertFrom} has not been converted to ${convertTo}`, 5000);
    };

    await testColumnAdded('smiles', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
    expect(DG.chem.isMolBlock(df.get('smiles_molblock', 0)));
    if (!DG.Test.isInBenchmark) {
      await testColumnAdded('smiles', DG.chem.Notation.Smiles, DG.chem.Notation.V3KMolBlock);
      expect(DG.chem.isMolBlock(df.get('smiles_v3Kmolblock', 0)));
      await testColumnAdded('smiles', DG.chem.Notation.Smiles, DG.chem.Notation.Smarts);
      expect(df.get('smiles_smarts', 0) === '[#6]-[#6](-[#6](=[#8])-[#8]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#7]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1)-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1');

      await testColumnAdded('smiles_molblock', DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
      expect(df.get('smiles_molblock_smiles', 0) === 'CC(C(=O)OCCCc1cccnc1)c1cccc(C(=O)c2ccccc2)c1');
      await testColumnAdded('smiles_molblock', DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock);
      expect(DG.chem.isMolBlock(df.get('smiles_molblock_v3Kmolblock', 0)));
      await testColumnAdded('smiles_molblock', DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
      expect(df.get('smiles_molblock_smarts', 0) === '[#6]-[#6](-[#6](=[#8])-[#8]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#7]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1)-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1');

      await testColumnAdded('smiles_v3Kmolblock', DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles);
      expect(df.get('smiles_v3Kmolblock_smiles', 0) === 'CC(C(=O)OCCCc1cccnc1)c1cccc(C(=O)c2ccccc2)c1');
      await testColumnAdded('smiles_v3Kmolblock', DG.chem.Notation.V3KMolBlock, DG.chem.Notation.MolBlock);
      expect(DG.chem.isMolBlock(df.get('smiles_v3Kmolblock_molblock', 0)));
      await testColumnAdded('smiles_v3Kmolblock', DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smarts);
      expect(df.get('smiles_v3Kmolblock_smarts', 0) === '[#6]-[#6](-[#6](=[#8])-[#8]-[#6]-[#6]-[#6]-[#6]1:[#6]:[#6]:[#6]:[#7]:[#6]:1)-[#6]1:[#6]:[#6]:[#6]:[#6](:[#6]:1)-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1');
    }
  }, {benchmark: true});
});
