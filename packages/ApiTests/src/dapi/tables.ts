import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, test} from '@datagrok-libraries/test/src/test';

category('Dapi: tables', () => {
  test('create and upload table', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('col', ['1', '2', '3'])]);
    df.name = 'dataframe';
    const id = await grok.dapi.tables.uploadDataFrame(df);
    await grok.dapi.tables.getTable(id);
  }, {stressTest: true,  owner: 'aparamonov@datagrok.ai'});
});
