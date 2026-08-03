import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

const TIMEOUT = 60 * 60 * 1000;

const BINDER_SEQUENCE = 'MKTAYIVKSHFSRQ';
const ASPIRIN_SMILES = 'CC(=O)OC1=CC=CC=C1C(=O)O';
const IBUPROFEN_SMILES = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O';
const CAFFEINE_SMILES = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C';
const PHENOL_SMILES = 'C1=CC=C(C=C1)O';
const TOLUENE_SMILES = 'CC1=CC=CC=C1';

category('Hosted API', () => {
  function expectFrame(result: any, rowCount: number, columns: string[]): DG.DataFrame {
    expect(result instanceof DG.DataFrame, true);
    const df = result as DG.DataFrame;
    expect(df.rowCount, rowCount);
    for (const name of columns)
      expect(df.columns.contains(name), true, `missing column "${name}"`);
    return df;
  }

  test('structureAndBinding', async () => {
    const table = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', [ASPIRIN_SMILES])]);
    const result = await grok.functions.call('Boltz1:BoltzStructureAndBinding',
      {table, ligands: table.getCol('smiles'), config: 'test-sab'});
    expectFrame(result, 1, ['binding_confidence']);
  }, {timeout: TIMEOUT});

  test('adme', async () => {
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('smiles', [ASPIRIN_SMILES, IBUPROFEN_SMILES, CAFFEINE_SMILES]),
    ]);
    const result = await grok.functions.call('Boltz1:BoltzAdme',
      {table, molecules: table.getCol('smiles')});
    expectFrame(result, 3, ['lipophilicity', 'permeability', 'solubility']);
  }, {timeout: TIMEOUT});

  test('smallMolecule.design', async () => {
    const result = await grok.functions.call('Boltz1:BoltzDesignSmallMolecules',
      {config: 'test-sm-design', numMolecules: 10});
    expect(result instanceof DG.DataFrame, true);
    expect((result as DG.DataFrame).rowCount > 0, true);
    expect((result as DG.DataFrame).columns.contains('smiles'), true);
  }, {timeout: TIMEOUT});

  test('smallMolecule.screen', async () => {
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('smiles', [ASPIRIN_SMILES, PHENOL_SMILES, TOLUENE_SMILES]),
    ]);
    const result = await grok.functions.call('Boltz1:BoltzScreenSmallMolecules',
      {table, molecules: table.getCol('smiles'), config: 'test-sm-screen'});
    expectFrame(result, table.rowCount, ['smiles']);
  }, {timeout: TIMEOUT});

  test('protein.design', async () => {
    const result = await grok.functions.call('Boltz1:BoltzDesignProteins',
      {config: 'test-protein-design', numProteins: 10});
    expect(result instanceof DG.DataFrame, true);
    expect((result as DG.DataFrame).rowCount > 0, true);
    expect((result as DG.DataFrame).columns.contains('sequence'), true);
  }, {timeout: TIMEOUT});

  test('protein.screen', async () => {
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('proteins', [BINDER_SEQUENCE, 'ACDEFGHIKLMNPQRSTVWY']),
    ]);
    const result = await grok.functions.call('Boltz1:BoltzScreenProteins',
      {table, proteins: table.getCol('proteins'), config: 'test-protein-screen'});
    expectFrame(result, table.rowCount, []);
  }, {timeout: TIMEOUT});
});
