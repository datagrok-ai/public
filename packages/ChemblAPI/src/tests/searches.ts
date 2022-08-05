import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, expect, delay } from '@datagrok-libraries/utils/src/test';
import { _package } from '../package-test';
import { chemblSimilaritySearch, chemblSubstructureSearch } from '../package';

const smiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';

category('Searches', () => {
  test('Similarity', async () => {
    const df = await grok.data.query(`${_package.name}:SimilaritySmileScore`, {'smile': smiles, 'score': 40});

    grok.shell.warning(df!.rowCount.toString());
    expect(df!.rowCount, 20);

    if (df != null)
      grok.shell.closeTable(df);
  });

  test('Substructure', async () => {
    const df = await grok.data.query(`${_package.name}:SubstructureSmile`, {'smile': smiles});
    grok.shell.warning(df!.rowCount.toString());
    expect(df!.rowCount, 20);

    if (df != null)
      grok.shell.closeTable(df);
  });
});
