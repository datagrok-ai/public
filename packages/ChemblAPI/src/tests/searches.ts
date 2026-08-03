import * as grok from 'datagrok-api/grok';

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {getData, SEARCH_TYPE} from '../utils';

const smiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';

category('Searches', () => {
  test('Similarity', async () => {
    const df = await getData(SEARCH_TYPE.SIMILARITY, smiles, 40);
    expect(df!.rowCount, 20);
    grok.shell.closeTable(df!);
  });

  test('Substructure', async () => {
    const df = await getData(SEARCH_TYPE.SUBSTRUCTURE, smiles);
    expect(df!.rowCount, 20);
    grok.shell.closeTable(df!);
  });
});
