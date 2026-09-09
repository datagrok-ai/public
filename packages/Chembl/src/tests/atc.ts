import * as grok from 'datagrok-api/grok';

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

category('ATC classification', () => {
  test('Four levels', async () => {
    const df = await grok.data.query(`${_package.name}:AtcClassification`, {
      'level1': 'ANTIINFECTIVES FOR SYSTEMIC USE', 'level2': 'ANTIMYCOTICS FOR SYSTEMIC USE',
      'level3': 'ANTIMYCOTICS FOR SYSTEMIC USE', 'level4': 'Triazole and tetrazole derivatives'});

    expect(df?.rowCount, 6);

    if (df != null)
      grok.shell.closeTable(df);
  });

  test('Level 1 and substructure', async () => {
    const df = await grok.data.query(`${_package.name}:AtcClassificationWithSubstructure`,
      {'level1': 'ANTIINFECTIVES FOR SYSTEMIC USE', 'substructure': 'Clc1ccccc1'});

    expect(df?.rowCount, 27);

    if (df != null)
      grok.shell.closeTable(df);
  });
});
