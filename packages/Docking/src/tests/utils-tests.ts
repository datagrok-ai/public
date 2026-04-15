import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {getFromPdb, buildComparisonTable} from '../utils/utils';

const SAMPLE_PDB = `REMARK  1  receptor. 1bdq J.
REMARK  2  binding energy.  -12.70
REMARK  3  intermolecular (1).  -15.68
REMARK  4  electrostatic.  -0.15
REMARK  5  ligand fixed.  -15.68
REMARK  6  ligand moving.  0.00
REMARK  7  total internal (2).  -1.55
REMARK  8  torsional free (3).  2.98
REMARK  9  unbound systems (4).  -1.55
ATOM      1  C1  LIG     1       0.000   0.000   0.000`;

const EXPECTED_VALUES: { [name: string]: number } = {
  'binding energy': -12.70, 'intermolecular (1)': -15.68,
  'electrostatic': -0.15, 'ligand fixed': -15.68,
  'ligand moving': 0.00, 'total internal (2)': -1.55,
  'torsional free (3)': 2.98, 'unbound systems (4)': -1.55,
};

category('Utils', () => {
  test('getFromPdb: parses energy values', async () => {
    const values = getFromPdb(SAMPLE_PDB);
    for (const [name, expected] of Object.entries(EXPECTED_VALUES))
      expect(values[name], expected);
  });

  test('getFromPdb: returns empty for no remarks', async () => {
    expect(Object.keys(getFromPdb('ATOM 1 C1 LIG 1 0 0 0')).length, 0);
  });

  test('buildComparisonTable: creates grid rows', async () => {
    const hovered = {'binding energy': -10.50, 'electrostatic': -2.36};
    const table = buildComparisonTable(EXPECTED_VALUES, hovered);
    expect(table instanceof HTMLElement, true);
    expect(table.querySelectorAll('[style*="grid"]').length > 0, true);
  });

  test('buildComparisonTable: handles missing hovered property', async () => {
    const table = buildComparisonTable(EXPECTED_VALUES, {'binding energy': -10.50});
    expect(table instanceof HTMLElement, true);
  });
}, {owner: 'oserhiienko@datagrok.ai'});
