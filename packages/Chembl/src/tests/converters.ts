import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, delay } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';

category('Converting queries', () => {
  test('Substructure pattern', async () => {
    const df = await grok.data.query(`${_package.name}:patternSubstructureSearch`, {'pattern': 'c1ccccc1', 'maxRows': 1000});

    expect(df!.rowCount, 1000);

    if (df != null)
      grok.shell.closeTable(df);
  });
});