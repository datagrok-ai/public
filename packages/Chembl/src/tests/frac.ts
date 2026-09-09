import * as grok from 'datagrok-api/grok';

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

category('FRAC classification', () => {
  test('By mechanism', async () => {
    const df = await grok.data.query(`${_package.name}:FracClassification`, {'mechanism': 'tubulin polymerization'});

    expect(df?.rowCount, 8);

    if (df != null)
      grok.shell.closeTable(df);
  });

  test('By mechanism and substructure', async () => {
    const df = await grok.data.query(`${_package.name}:FracClassificationWithSubstructure`,
      {'mechanism': 'C14-demethylase in sterol biosynthesis (erg11/cyp51)', 'substructure': 'Clc1ccccc1'});

    expect(df?.rowCount, 26);

    if (df != null)
      grok.shell.closeTable(df);
  });
});
