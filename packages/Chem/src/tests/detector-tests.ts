import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {dfFromColWithOneCategory, loadFileAsBytes, loadFileAsText} from './utils';
import {_importSdf} from '../open-chem/sdf-importer';

category('detectors', () => {
  test('detectSmiles', async () => {
    const df = DG.DataFrame.fromCsv(await loadFileAsText('sar-small.csv'));
    let col = df.columns.byName('smiles');
    await grok.data.detectSemanticTypes(df);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.tags[DG.TAGS.UNITS], DG.UNITS.Molecule.SMILES);

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
    expect(col.tags[DG.TAGS.UNITS], DG.UNITS.Molecule.SMILES);

    const df3 = DG.DataFrame.fromCsv(testSmiles2);
    col = df3.columns.byName('smiles');
    col.set(2, 'not a molstring');
    col.set(4, 'not a molstring');
    await grok.data.detectSemanticTypes(df3);
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.tags[DG.TAGS.UNITS], DG.UNITS.Molecule.SMILES);
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
    expect(col.tags[DG.TAGS.UNITS], DG.UNITS.Molecule.SMILES);

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
    expect(col.tags[DG.TAGS.UNITS], DG.UNITS.Molecule.MOLBLOCK);
  });

  test('detectMolblockSDF', async () => {
    const [df] = _importSdf(await loadFileAsBytes('mol1K.sdf'));
    await grok.data.detectSemanticTypes(df);
    const col = df.columns.byName('molecule');
    expect(col.semType, DG.SEMTYPE.MOLECULE);
    expect(col.tags[DG.TAGS.UNITS], DG.UNITS.Molecule.MOLBLOCK);
  });

  test('spgi', async () => {
    const df = DG.DataFrame.fromCsv(await loadFileAsText('tests/spgi-100.csv'));
    await grok.data.detectSemanticTypes(df);
    expect(df.col('structure')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('primary scaffold')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('core')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r1')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r2')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r3')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r100')!.semType, DG.SEMTYPE.MOLECULE);
    expect(df.col('r101')!.semType, DG.SEMTYPE.MOLECULE);
  });
});
