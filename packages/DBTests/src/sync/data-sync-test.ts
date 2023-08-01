import {awaitCheck, category, expect, test} from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('data sync', () => {
  test('grok connect streaming and data sync', async () => {
    grok.functions.eval('DbTests:PostgresqlTestDataSync()').then(() => {});
    await awaitCheck(() => grok.shell.tv?.table?.name == 'result', 'Query first batch timeout', 30000);
    expect(grok.shell.tv?.table?.tags['.script'], 'Result = Dbtests:PostgresqlTestDataSync()');
    grok.shell.closeTable(grok.shell.tv!.table!);
    grok.shell.tv?.close();
  });
});
