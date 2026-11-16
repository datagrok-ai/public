import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {dfFromColWithOneCategory, loadFileAsBytes, loadFileAsText} from './utils';
import {_importSdf} from '../open-chem/sdf-importer';

category('detectors', () => {
  test('detectSmiles', async () => {
    const df = DG.DataFrame.fromCsv(await loadFileAsText('sar-small.csv'));
    let col = df.columns.byName('smiles');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.meta.units, DG.UNITS.Molecule.SMILES);

    const testSmiles2 = `smiles
O=C1CN=C(c2ccccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`;
    const df2 = DG.DataFrame.fromCsv(testSmiles2);
    col = df2.columns.byName('smiles');
    await grok.data.detectSemanticTypes(df2);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.meta.units, DG.UNITS.Molecule.SMILES);

    const df3 = DG.DataFrame.fromCsv(testSmiles2);
    col = df3.columns.byName('smiles');
    col.set(2, 'not a molstring');
    col.set(4, 'not a molstring');
    await grok.data.detectSemanticTypes(df3);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.meta.units, DG.UNITS.Molecule.SMILES);
  });

  test('detectShortSmiles', async () => {
    let df = await dfFromColWithOneCategory('result', 'NO', 10);
    let col = df.columns.byName('result');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType != DG.SEMTYPE.MOLECULE, true);

    df = await dfFromColWithOneCategory('smiles', 'NO', 10);
    col = df.columns.byName('smiles');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.meta.units, DG.UNITS.Molecule.SMILES);

    df = await dfFromColWithOneCategory('smiles', 'OK', 10);
    col = df.columns.byName('smiles');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType != DG.SEMTYPE.MOLECULE, true);
  });

  test('detectMolblock', async () => {
    const df = DG.DataFrame.fromCsv(await loadFileAsText('molblocks.csv'));
    const col = df.columns.byName('scaffold');
    col.name = 'not familiar name';
    await grok.data.detectSemanticTypes(df);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.meta.units, DG.UNITS.Molecule.MOLBLOCK);
  });

  test('detectMolblockSDF', async () => {
    const [df] = _importSdf(await loadFileAsBytes('mol1K.sdf'));
    await grok.data.detectSemanticTypes(df);
    const col = df.columns.byName('molecule');
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.meta.units, DG.UNITS.Molecule.MOLBLOCK);
  });

  test('detectReactionTemplates', async () => {
    // Chem's own enumerator templates hold SMARTS with `!`, `&` and `~`: a reaction column, not text
    const templates = DG.DataFrame.fromCsv(await loadFileAsText('enumerations/reactions.csv'));
    await grok.data.detectSemanticTypes(templates);
    expect(templates.col('reaction_smarts')!.semType, DG.SEMTYPE.CHEMICAL_REACTION);

    // eight reactions in ten values are a reaction column: the detector asks for 80 %
    const steps = DG.DataFrame.fromColumns([DG.Column.fromStrings('steps', [
      'CO>>COC', 'CC(=O)O.N>>CC(=O)N', 'CBr.c1ccccc1B(O)O>>Cc1ccccc1', 'CCO>>CC=O', 'CC=O>>CC(=O)O',
      'c1ccccc1>>Brc1ccccc1', 'CCBr>>CC[Mg]Br', 'C=C>>CC', 'not recorded', 'see the notebook'])]);
    await grok.data.detectSemanticTypes(steps);
    expect(steps.col('steps')!.semType, DG.SEMTYPE.CHEMICAL_REACTION);

    // plain text with `>` in it (and `!`, `&`, `~`) stays text
    const notes = DG.DataFrame.fromColumns([DG.Column.fromStrings('notes', [
      'yield > 90 %', 'R&D > operations', 'heat to 80 C >> cool', 'ratio ~ 2:1 > 1:1', 'done!>next',
      'see A > B > C', 'pH>7', 'step 1 >> step 2', 'warm & stir > 2 h', 'T > Tm!'])]);
    await grok.data.detectSemanticTypes(notes);
    expect(notes.col('notes')!.semType !== DG.SEMTYPE.CHEMICAL_REACTION, true);
  });

  test('detectReactionUnits', async () => {
    // a reaction column's units name the notation of its values: smiles, smarts, rxn or route
    const unitsOf = async (name: string, values: string[]) => {
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(name, values)]);
      await grok.data.detectSemanticTypes(df);
      expect(df.col(name)!.semType, DG.SEMTYPE.CHEMICAL_REACTION);
      return df.col(name)!.meta.units;
    };
    expect(await unitsOf('reaction', ['CCO>>CC=O', 'CC(=O)O.OCC>OS(=O)(=O)O>CC(=O)OCC.O',
      'CC>>CO |(0,0,;1,0,;0,0,;1,0,)|']), 'smiles');
    expect(await unitsOf('reaction', ['CCO>>CC=O', '[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]',
      '[c:1][Br,I].[c:2]B(O)O>>[c:1][c:2]']), 'smarts');
    const rxn = await loadFileAsText('tests/reactions-crux.csv');
    const blocks = DG.DataFrame.fromCsv(rxn).col('reaction')!.toList().filter((s: string) => s.startsWith('$RXN'));
    expect(await unitsOf('reaction', blocks), 'rxn');
    expect(await unitsOf('reaction', ['CCO>>CC=O', 'CCO>>CC=O>>CC(=O)O', 'CO>>COC--**--COC.Cl>>CCl']), 'route');
    // routes whose steps carry agents: the detector judges a route by its first step
    expect(await unitsOf('reaction', ['CCO>O=[Cr](=O)=O>CC=O--**--CC=O.N>>CC=N',
      'Oc1ccccc1.O=C=O>[Na+].[OH-]>O=C(O)c1ccccc1O--**--' +
      'O=C(O)c1ccccc1O.CC(=O)OC(C)=O>O=P(O)(O)O>CC(=O)Oc1ccccc1C(=O)O']),
    'route');
    const templates = DG.DataFrame.fromCsv(await loadFileAsText('enumerations/reactions.csv'));
    await grok.data.detectSemanticTypes(templates);
    expect(templates.col('reaction_smarts')!.meta.units, 'smarts');
    const samples = DG.DataFrame.fromCsv(await loadFileAsText('test-reactions.csv'));
    await grok.data.detectSemanticTypes(samples);
    expect(samples.col('reaction')!.meta.units, 'route');
  });

  test('spgi', async () => {
    const df = DG.DataFrame.fromCsv(await loadFileAsText('tests/spgi-100.csv'));
    await grok.data.detectSemanticTypes(df);
    //'Primary scaffold' column has been removed from spgi version 2
    expect(df.col('Structure')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('core')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r1')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r2')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r3')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r100')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r101')!.semType, DG.SEMTYPE.MOLECULE);
  });
});
