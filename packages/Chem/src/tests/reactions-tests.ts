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
import {ChemicalEngine} from '../engine/chem-engine';
import {expectFallback, inBothEngines, inEngine} from './engine-utils';
import {cruxLoadFailed, resetCruxLoaderForTests} from '../engine/crux-module';
import {lastReactionEngine, setCruxUnreadableReactantsForTests} from '../utils/reactions/reaction-engine';
import {validateReactionSmarts} from '../utils/reactions/reaction-storage';
import {openReactionEditor} from '../utils/reactions/reaction-editor';

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

  test('Crux mode: the products of RDKit mode, from Crux', async () => {
    // the transformation, two-component and deprotection reactions run on the Crux module in Crux mode
    // (utils/reactions/reaction-engine.ts), molecule by molecule, with RDKit behind it: the same products.
    // lastReactionEngine tells the engine the reaction was set up on. The deprotection's `C1CC(` is read by neither
    // engine: it fails alike in both modes, no fallback (it was counted as one, RDKit failing it too).
    const molecules = (await readDataframe('smiles.csv')).col('canonical_smiles')!.toList().slice(0, 300);
    const acids = ['OC(=O)c1ccccc1', 'CC(=O)O', 'OC(=O)C1CCCCC1', 'OC(=O)c1cccnc1', 'Cc1ccccc1'];
    const amines = ['Nc1ccccc1', 'NC1CCCCC1', 'CCN', 'COc1ccc(N)cc1', 'COc1ccccc1'];
    const fmoc = DG.Column.fromStrings('smiles', ['O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21',
      'O=C(O)C1CCCN1C(=O)OCC1c2ccccc2-c2ccccc21', 'CC(C)C(N)C(=O)O', 'C1CC(']);
    const run = async (engine: ChemicalEngine) => {
      const hydrolysis = await runTransformationReaction(rdkit, smartsOf('Ester Hydrolysis (Saponification)'),
        molecules);
      expect(lastReactionEngine(), engine, 'the transformation\'s engine');
      const coupling = await runTwoComponentReaction(rdkit, smartsOf('Amide Coupling'), acids, amines,
        {mode: 'matrix'});
      expect(lastReactionEngine(), engine, 'the two-component reaction\'s engine');
      const deprotected = (await addDeprotectedColumn(DG.DataFrame.fromColumns([fmoc.clone()]),
        fmoc, defaultFragment, rdkit))!.toList();
      expect(lastReactionEngine(), engine, 'the deprotection\'s engine');
      return {hydrolysis: canonical(hydrolysis.products), errors: hydrolysis.errors,
        coupling: canonical(coupling.products), deprotected: canonical(deprotected.slice(0, 3)).concat(deprotected[3])};
    };
    const {rdkit: rdkitMode, crux: cruxMode} = await inBothEngines(run,
      {noFallback: ['transformation reaction', 'two-component reaction', 'deprotection']});
    expect(cruxMode.hydrolysis.filter((p, i) => p !== rdkitMode.hydrolysis[i]).length, 0, 'hydrolysis products');
    expect(cruxMode.errors, rdkitMode.errors);
    expectArray(cruxMode.coupling, rdkitMode.coupling);
    expect(cruxMode.coupling.filter((p) => p).length, 16, 'amide couplings');
    expectArray(cruxMode.deprotected, rdkitMode.deprotected);
    expectArray(cruxMode.deprotected, canonical(['NCC(=O)O', 'O=C(O)C1CCCN1', 'CC(C)C(N)C(=O)O']).concat('C1CC('));
    // a molecule Crux cannot read gets RDKit's product, the fallback counted
    const esters = ['CCOC(C)=O', 'CC(=O)OC'];
    const hydrolysis = async () => canonical((await runTransformationReaction(rdkit,
      smartsOf('Ester Hydrolysis (Saponification)'), esters)).products);
    const rdkitProducts = await inEngine(ChemicalEngine.RDKit, hydrolysis);
    await inEngine(ChemicalEngine.Crux, async () => {
      setCruxUnreadableReactantsForTests([esters[0]]);
      try {
        const {result, reason} = await expectFallback('transformation reaction', hydrolysis);
        expectArray(result, rdkitProducts);
        // Crux's own words for it, not a type error out of its glue ("expected instance of w": a null molecule
        // handed on)
        expect(reason, `no molecule read from ${esters[0]}`, 'why RDKit did the molecule');
      } finally {
        setCruxUnreadableReactantsForTests(null);
      }
    });
  });

  test('Crux mode: a molecule only RDKit reads goes to RDKit, named; one neither reads fails alike, no fallback',
    async () => {
      const amide = smartsOf('Amide Coupling');
      const acids = ['OC(=O)c1ccccc1', 'CC(=O)O'];
      const amines = ['Nc1ccccc1', 'CCN'];
      const couple = async () => canonical((await runTwoComponentReaction(rdkit, amide, acids, amines,
        {mode: 'pairwise'})).products);
      const rdkitProducts = await inEngine(ChemicalEngine.RDKit, couple);
      // the two-component reaction: a reactant Crux reads none from is RDKit's, the fallback naming it
      await inEngine(ChemicalEngine.Crux, async () => {
        setCruxUnreadableReactantsForTests([amines[1]]);
        try {
          const {result, reason} = await expectFallback('two-component reaction', couple);
          expectArray(result, rdkitProducts);
          expect(reason, `no molecule read from ${amines[1]}`, 'why RDKit did the pair');
        } finally {
          setCruxUnreadableReactantsForTests(null);
        }
      });
      // a molecule neither engine reads: no product and an error in both modes, and no fallback in Crux mode (RDKit
      // failing it too was counted as one)
      const hydrolysis = smartsOf('Ester Hydrolysis (Saponification)');
      const run = async () => {
        const one = await runTransformationReaction(rdkit, hydrolysis, ['CCOC(C)=O', 'C1CC(']);
        const two = await runTwoComponentReaction(rdkit, amide, ['CC(=O)O', 'C1CC('], ['CCN', 'CCN'],
          {mode: 'pairwise'});
        return {hydrolysed: canonical(one.products.slice(0, 1)), two: canonical(two.products), errors: two.errors};
      };
      const {rdkit: rdkitMode, crux: cruxMode} = await inBothEngines(run,
        {noFallback: ['transformation reaction', 'two-component reaction']});
      expect(JSON.stringify(cruxMode), JSON.stringify(rdkitMode), 'the same products and errors');
      expect(cruxMode.errors, 1, 'the two-component reaction\'s error: the molecule neither engine reads');
    });

  test('break it: a product RDKit 2026 cannot write leaves Crux running', async () => {
    // RDKit 2026.03 throws writing the bromination product of a trans-1,4-cyclohexane ("no other ring atoms found");
    // Crux trapped there, which sent the rest of the session to RDKit (the 100K comparison, 2026-10-23): now it throws,
    // and RDKit writes that molecule's product
    const molecules = ['O[C@@H]1CC[C@H](c2ccccc2)CC1', 'c1ccncc1'];
    const brominate = async () => (await runTransformationReaction(rdkit, '[cH:1]>>[c:1]Br', molecules)).products;
    try {
      const {rdkit: rdkitProducts, crux: cruxProducts} = await inBothEngines(brominate);
      expect(cruxLoadFailed(), false, 'the Crux module did not trap');
      expectArray(cruxProducts, rdkitProducts);
    } finally {
      if (cruxLoadFailed())
        resetCruxLoaderForTests(null);
    }
  });

  test('reaction editor: a reaction with agents opens and is accepted', async () => {
    // reactants>agents>products, as the grid draws it; the editor refused it ("must contain exactly one >>")
    const esterification = 'CC(=O)O.OCC>OS(=O)(=O)O>CC(=O)OCC.O';
    const check = validateReactionSmarts(rdkit, esterification);
    expect(check.valid, true, check.error);
    expect(check.numReactants, 2);
    expect(check.numProducts, 2);
    expect(validateReactionSmarts(rdkit, 'CCO>>CC=O>>CC(=O)O').valid, false, 'a route is not one named reaction');

    await openReactionEditor(rdkit, {preset: {id: 'test-agents', name: 'Fischer esterification',
      reactionSmarts: esterification, mode: 'two-component', category: 'Test'}});
    const dialog = DG.Dialog.getOpenDialogs().find((d) => d.title === 'Edit Reaction')!;
    try {
      expect(!!dialog, true, 'the editor opened');
      expect(dialog.root.textContent!.includes('Valid  |  Reactants: 2  |  Products: 2'), true,
        dialog.root.textContent!);
      expect(dialog.getButton('OK').classList.contains('d4-disabled'), false, 'SAVE is enabled');
    } finally {
      dialog?.close();
    }
  });

  test('deprotect: the default fragment takes Fmoc off', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', [
      'O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21', 'CC(NC(=O)OCC1c2ccccc2-c2ccccc21)C(=O)O',
      'O=C(O)C1CCCN1C(=O)OCC1c2ccccc2-c2ccccc21', 'CC(C)C(N)C(=O)O'])]);
    const deprotected = (await addDeprotectedColumn(df, df.col('smiles')!, defaultFragment, rdkit))!;
    expectArray(canonical(deprotected.toList()),
      canonical(['NCC(=O)O', 'CC(N)C(=O)O', 'O=C(O)C1CCCN1', 'CC(C)C(N)C(=O)O']));
  });
});
