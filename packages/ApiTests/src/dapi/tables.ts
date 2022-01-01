import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('Dapi: tables', () => {

  test('Dapi: tables - create and upload table', async () => {
    let df = DG.DataFrame.fromColumns([DG.Column.fromStrings('col', ['1', '2', '3'])]);
    df.name = 'dataframe';
    let id = await grok.dapi.tables.uploadDataFrame(df);
    await grok.dapi.tables.getTable(id);
  });

});
