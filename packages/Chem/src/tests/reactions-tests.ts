import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {PackageFunctions} from '../package';
import {removeWaterAndSalts, runTransformationReaction, runTwoComponentReaction} from '../utils/reactions/reactions';
import {NamedReaction} from '../utils/reactions/types';
import {addDeprotectedColumn, defaultFragment} from '../analysis/deprotect';
import {loadFileAsText, readDataframe} from './utils';

category('Reactions', () => {
  let rdkit: RDModule;
  let reactions: NamedReaction[];

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    rdkit = PackageFunctions.getRdKitModule();
    reactions = JSON.parse(await loadFileAsText('reactions/default-reactions.json'));
  });

  const smartsOf = (name: string) => reactions.find((r) => r.name === name)!.reactionSmarts;

  const canonical = (molecules: string[]) => molecules.map((s) => {
    if (!s)
      return '';
    const mol = rdkit.get_mol(s);
    try {
      return mol.get_smiles();
    } finally {
      mol.delete();
    }
  });

  const changedRows = (a: string[], b: string[]) => a.map((_, i) => i).filter((i) => a[i] !== b[i]);

  test('remove water and salts', async () => {
    const df = await readDataframe('smiles.csv');
    await grok.data.detectSemanticTypes(df);
    const original = canonical(df.col('canonical_smiles')!.toList());
    await PackageFunctions.removeWaterAndSaltsTopMenu(df, df.col('canonical_smiles')!);
    const desalted = df.col('Desalted(canonical_smiles)')!;
    expect(desalted.semType, DG.SEMTYPE.MOLECULE);
    expect(desalted.stats.missingValueCount, 0);
    const after = canonical(desalted.toList());
    expectArray(changedRows(after, original), [938, 939, 948, 949, 950, 951, 952, 953, 954, 955, 956, 957, 959]);
    for (let i = 0; i < after.length; i++) {
      expect(after[i].split('.').length <= original[i].split('.').length, true, `row ${i} gained a fragment`);
      for (const fragment of after[i].split('.'))
        expect(['O', 'Cl', 'Br', 'I', 'F'].includes(fragment), false, `row ${i} keeps ${fragment}`);
    }
    const halideAnions = after.map((_, i) => i).filter((i) => /\[(Cl|Br|I)-\]/.test(after[i]));
    expectArray(halideAnions, [51, 937, 940, 941, 942, 943, 944, 945, 947, 958]);

    const columns = df.columns.length;
    await PackageFunctions.removeWaterAndSaltsTopMenu(df, df.col('canonical_smiles')!);
    expect(df.columns.length, columns);
  });

  test('transformation: ester hydrolysis', async () => {
    const df = await readDataframe('smiles.csv');
    const molecules = df.col('canonical_smiles')!.toList();
    const res = await runTransformationReaction(rdkit, smartsOf('Ester Hydrolysis (Saponification)'), molecules);
    const desalted = canonical(await removeWaterAndSalts(molecules));
    const products = canonical(res.products);
    expect(changedRows(products, desalted).length, 244);
    expect(products[0], 'CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1');
    expect(res.errors, 0);
  });

  test('two-component: amide coupling', async () => {
    const acids = ['OC(=O)c1ccccc1', 'CC(=O)O', 'CCC(=O)O', 'CCCC(=O)O', 'OC(=O)C1CCCCC1', 'Cc1ccc(cc1)C(=O)O',
      'OC(=O)c1ccc(Cl)cc1', 'OC(=O)Cc1ccccc1', 'OC(=O)c1cccnc1', 'OC(=O)c1ccco1', 'Cc1ccccc1', 'OC(=O)c1ccccc1'];
    const amines = ['Nc1ccccc1', 'NCc1ccccc1', 'NC1CCCCC1', 'CCCCN', 'Cc1ccc(N)cc1', 'Nc1ccc(Cl)cc1', 'NCCc1ccccc1',
      'Nc1cccnc1', 'CCN', 'COc1ccc(N)cc1', 'Nc1ccccc1', 'COc1ccccc1'];
    const res = await runTwoComponentReaction(rdkit, smartsOf('Amide Coupling'), acids, amines, {mode: 'pairwise'});
    expectArray(canonical(res.products), [
      'O=C(Nc1ccccc1)c1ccccc1', 'CC(=O)NCc1ccccc1', 'CCC(=O)NC1CCCCC1', 'CCCCNC(=O)CCC',
      'Cc1ccc(NC(=O)C2CCCCC2)cc1', 'Cc1ccc(C(=O)Nc2ccc(Cl)cc2)cc1', 'O=C(NCCc1ccccc1)c1ccc(Cl)cc1',
      'O=C(Cc1ccccc1)Nc1cccnc1', 'CCNC(=O)c1cccnc1', 'COc1ccc(NC(=O)c2ccco2)cc1', '', '']);
    expect(res.noProducts, 2);
  });

  test('deprotect: the default fragment takes Fmoc off', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', [
      'O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21', 'CC(NC(=O)OCC1c2ccccc2-c2ccccc21)C(=O)O',
      'O=C(O)C1CCCN1C(=O)OCC1c2ccccc2-c2ccccc21', 'CC(C)C(N)C(=O)O'])]);
    const deprotected = addDeprotectedColumn(df, df.col('smiles')!, defaultFragment, rdkit)!;
    expectArray(canonical(deprotected.toList()),
      canonical(['NCC(=O)O', 'CC(N)C(=O)O', 'O=C(O)C1CCCN1', 'CC(C)C(N)C(=O)O']));
  });
});
