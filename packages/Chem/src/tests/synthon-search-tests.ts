import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {readDataframe} from './utils';

const MOLECULE = 'Cc1ccc2c(C(NC)=O)n[nH]c2c1';
const SPACE = 'Syntons_5567.csv';

function compareDfColumns(result: DG.DataFrame, expected: DG.DataFrame): void {
  for (const col of expected.columns)
    expect(result.columns.contains(col.name), true, `missing column: ${col.name}`);
}

function compareProductColumn(result: DG.DataFrame, expected: DG.DataFrame): void {
  const resProducts = new Set<string>();
  const resCol = result.getCol('product');
  for (let i = 0; i < result.rowCount; i++)
    resProducts.add(resCol.get(i));

  const expCol = expected.getCol('product');
  for (let i = 0; i < expected.rowCount; i++)
    expect(resProducts.has(expCol.get(i)), true, `missing product: ${expCol.get(i)}`);
}

category('synthon search', () => {
  test('substructure search', async () => {
    const result: DG.DataFrame = await grok.functions.call('Chem:synthonSearchFunc', {
      spaceName: SPACE,
      molecule: MOLECULE,
      maxHits: 100,
      searchType: 'substructure',
    });
    expect(result instanceof DG.DataFrame, true, 'result should be a DataFrame');
    expect(result.rowCount > 0, true, 'should return at least one hit');

    const expected = await readDataframe('tests/Synthon_Substructure_Search_Results.csv');
    compareDfColumns(result, expected);
    expect(result.rowCount, expected.rowCount, 'row count mismatch');
    compareProductColumn(result, expected);
  }, {timeout: 120000, skipReason: 'GROK-19863'});

  test('similarity search', async () => {
    const result: DG.DataFrame = await grok.functions.call('Chem:synthonSearchFunc', {
      spaceName: SPACE,
      molecule: MOLECULE,
      maxHits: 100,
      searchType: 'similarity',
      similarityCutoff: 0.5,
    });
    expect(result instanceof DG.DataFrame, true, 'result should be a DataFrame');
    expect(result.rowCount > 0, true, 'should return at least one hit');

    const expected = await readDataframe('tests/Synthon_Similarity_Search_Results.csv');
    compareDfColumns(result, expected);
    expect(result.rowCount, expected.rowCount, 'row count mismatch');
    compareProductColumn(result, expected);

    const simCol = result.getCol('similarity');
    for (let i = 0; i < result.rowCount; i++)
      expect(simCol.get(i) >= 0.5, true, `similarity below cutoff at row ${i}`);
  }, {timeout: 120000, skipReason: 'GROK-19863'});
});
