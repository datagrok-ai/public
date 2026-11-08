import {before, category, test, expect} from '@datagrok-libraries/test/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {readDataframe} from './utils';

category('save as csv', async () => {
  let df: DG.DataFrame;

  before(async () => {
    df = await readDataframe('tests/spgi-100.csv');
    await grok.data.detectSemanticTypes(df);
  });

  test('molblocks exported as smiles read back as molecules', async () => {
    const structure = df.col('Structure')!;
    expect(structure.meta.units, DG.UNITS.Molecule.MOLBLOCK);
    const csv = await df.toCsvEx({moleculesAsSmiles: true});
    expect(csv.includes('M  END'), false);

    const back = DG.DataFrame.fromCsv(csv);
    await grok.data.detectSemanticTypes(back);
    const smiles = back.col('Structure')!;
    expect(back.rowCount, df.rowCount);
    expect(smiles.semType, DG.SEMTYPE.MOLECULE);
    expect(smiles.meta.units, DG.UNITS.Molecule.SMILES);
    expect(smiles.stats.missingValueCount, 0);
    expect(smiles.categories.length, structure.categories.length);
  });

  test('molblocks stay molblocks without the option', async () => {
    const csv = await df.toCsvEx({moleculesAsSmiles: false});
    expect(csv.includes('M  END'), true);
  });

  test('an empty molecule exports as an empty field', async () => {
    const structures = df.col('Structure')!;
    const t = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', ['1', '2', '3']),
      DG.Column.fromStrings('structure', [structures.get(0), '', structures.get(1)]),
      DG.Column.fromStrings('tag', ['a', 'b', 'c']),
    ]);
    await grok.data.detectSemanticTypes(t);
    const csv = await t.toCsvEx({moleculesAsSmiles: true});
    expect(csv.includes('M  END'), false);
    expect(csv.includes('\n2,,b'), true);
    expect(csv.includes('null'), false);
  });
});
