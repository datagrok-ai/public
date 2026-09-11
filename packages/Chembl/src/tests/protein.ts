import * as grok from 'datagrok-api/grok';

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';

category('Protein activity', () => {
  test('Targets containing protein', async () => {
    const df = await grok.data.query(`${_package.name}:CompoundActivityDetailsForAllTargetsContainingProtein`,
      {'protein': 'P08172'});

    expect(df?.rowCount, 15334);

    if (df != null)
      grok.shell.closeTable(df);
  }, {timeout: 60000});
});
