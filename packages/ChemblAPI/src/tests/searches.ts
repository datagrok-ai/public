import * as grok from 'datagrok-api/grok';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import { getData, SEARCH_TYPE } from '../package';

const smiles = 'O=C1CN=C(c2ccccc2N1)C3CCCCC3';

category('Searches', () => {
  test('Similarity', async () => {
    const df = getData(SEARCH_TYPE.SIMILARITY, smiles, 40);
    expect(df!.rowCount, 19);

    if (df != null)
      grok.shell.closeTable(df);
  });

  test('Substructure', async () => {
    const df = getData(SEARCH_TYPE.SUBSTRUCTURE, smiles, 40);
    expect(df!.rowCount, 19);

    if (df != null)
      grok.shell.closeTable(df);
  });
});
