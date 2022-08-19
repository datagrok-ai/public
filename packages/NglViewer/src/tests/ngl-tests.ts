import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, delay } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';

category('Cartridge usage', () => {
  // test('Similarity pattern', async () => {
  //   const df = await grok.data.query(`${_package.name}:patternSimilaritySearch`, {'pattern': 'c1ccc(O)cc1', 'maxRows': 1000});

  //   expect(df!.rowCount, 25);

  //   if (df != null)
  //     grok.shell.closeTable(df);
  // });
});
