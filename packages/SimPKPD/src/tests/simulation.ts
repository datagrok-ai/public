import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, delay } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';
import { simulatePKPD } from '../package';

category('Simulation', () => {
  test('Simulation', async () => {
    
    const df = await simulatePKPD(10000.0, 10, 12, 0.3, 2.0, 4.0, 1.0, 30.0, 0.2, 8.0);

    expect(df!.rowCount, 1010);    
  });
});

