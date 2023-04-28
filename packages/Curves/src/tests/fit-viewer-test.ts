import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';


category('Fit', () => {
  test('test', async () => {
    expect('test', 'test');
  });
});
