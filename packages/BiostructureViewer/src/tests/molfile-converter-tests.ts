import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {molfileV3KToV2K} from '../utils/molfile-converter';

const benzeneV3K = `v3000 benzene
  Datagrok

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 1.5 0.0 0
M  V30 2 C 1.299 0.75 0.0 0
M  V30 3 C 1.299 -0.75 0.0 0
M  V30 4 C 0.0 -1.5 0.0 0
M  V30 5 C -1.299 -0.75 0.0 0
M  V30 6 C -1.299 0.75 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 END BOND
M  V30 END CTAB
M  END
`;

category('molfileConverter', () => {
  test('molfileV3KToV2K', async () => {
    const v2k = molfileV3KToV2K(benzeneV3K);
    const lines = v2k.split('\n');
    expect(lines[3], '  6  6  0  0  0  0  0  0  0  0999 V2000');
    expect(lines[4], '    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0');
    expect(lines[8], '   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0');
    expect(lines[10], '  1  2  2  0');
    expect(lines[15], '  6  1  1  0');
    expect(lines[16], 'M  END');

    expect(MolfileHandler.isMolfileV2K(v2k), true);
    const check = MolfileHandler.getInstance(v2k);
    expect(check.atomCount, 6);
    expect(check.bondCount, 6);
    expect(check.atomTypes.every((t) => t === 'C'), true);
  });
});
