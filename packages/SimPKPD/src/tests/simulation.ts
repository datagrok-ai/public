import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, delay } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';
import { simulate } from '../package';

category('Simulation', () => {
  test('Simulation', async () => {

    const df = await simulate(1000, 12, '2 compartment PK', 2, 0.3, 4, 30, 1, 0.2, 8);

    expect(df!.rowCount, 200);

    if (df != null)
      grok.shell.closeTable(df);
  });
});
