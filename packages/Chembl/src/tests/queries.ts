import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, delay } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';

const smiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';

category('Searches', () => {
  test('Similarity pattern', async () => {
    const df = await grok.data.query(`${_package.name}:patternSimilaritySearch`, {'pattern': 'c1ccc(O)cc1', 'maxRows': 1000});

    expect(df!.rowCount, 25);

    if (df != null)
      grok.shell.closeTable(df);
  });
});
