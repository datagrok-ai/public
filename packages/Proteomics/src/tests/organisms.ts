import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {resolveOrganismCode, detectOrganismCode, ORGANISM_LIST} from '../utils/organisms';

function dfWithOrganismColumn(name: string, values: string[]): DG.DataFrame {
  return DG.DataFrame.fromColumns([DG.Column.fromStrings(name, values)]);
}

category('Organisms', () => {
  test('resolveOrganismCode maps scientific names (incl. strain-qualified) to codes', async () => {
    expect(resolveOrganismCode('Homo sapiens'), 'hsapiens');
    expect(resolveOrganismCode('Rattus norvegicus'), 'rnorvegicus');
    // Strain-qualified values still resolve via the scientific-name prefix.
    expect(resolveOrganismCode('Escherichia coli (strain K12)'), 'ecoli');
    expect(resolveOrganismCode('Saccharomyces cerevisiae (strain ATCC 204508 / S288c)'), 'scerevisiae');
  });

  test('resolveOrganismCode returns undefined for unsupported / empty', async () => {
    expect(resolveOrganismCode('Sus scrofa') === undefined, true);
    expect(resolveOrganismCode('') === undefined, true);
    expect(resolveOrganismCode(null) === undefined, true);
    expect(resolveOrganismCode(undefined) === undefined, true);
  });

  test('every ORGANISM_LIST display resolves back to its own code', async () => {
    for (const o of ORGANISM_LIST)
      expect(resolveOrganismCode(o.display), o.code);
  });

  test('detectOrganismCode: single-species PG.Organisms column → that code', async () => {
    const df = dfWithOrganismColumn('PG.Organisms',
      ['Rattus norvegicus', 'Rattus norvegicus', 'Rattus norvegicus']);
    expect(detectOrganismCode(df), 'rnorvegicus');
  });

  test('detectOrganismCode: works on a generically named "Organism" column', async () => {
    const df = dfWithOrganismColumn('Organism', ['Homo sapiens', 'Homo sapiens']);
    expect(detectOrganismCode(df), 'hsapiens');
  });

  test('detectOrganismCode: mixed empties/unknowns but one species → that code', async () => {
    const df = dfWithOrganismColumn('PG.Organisms', ['', 'Homo sapiens', 'Sus scrofa', '']);
    expect(detectOrganismCode(df), 'hsapiens');
  });

  test('detectOrganismCode: multi-species (HYE-style mix) → undefined (do not guess)', async () => {
    const df = dfWithOrganismColumn('PG.Organisms',
      ['Homo sapiens', 'Escherichia coli (strain K12)', 'Saccharomyces cerevisiae']);
    expect(detectOrganismCode(df) === undefined, true);
  });

  test('detectOrganismCode: no organism column → undefined', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('Gene', ['CDK1', 'CCNB1'])]);
    expect(detectOrganismCode(df) === undefined, true);
  });
});
