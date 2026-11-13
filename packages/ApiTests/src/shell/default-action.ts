import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

category('ObjectHandler: runDefaultAction', () => {
  test('a CSV file opens as a table, as a double-click in Browse does', async () => {
    const [file] = (await grok.dapi.files.list('System:DemoFiles/', false, 'csv'))
      .filter((f) => f.fileName === 'demog.csv');
    expect(file != null, true, 'System:DemoFiles/demog.csv is there');
    const before = new Set(Array.from(grok.shell.tableViews).map((v) => v.id));
    let opened: DG.TableView | undefined;
    try {
      DG.ObjectHandler.runDefaultAction(file);
      await waitUntil(() => (opened = Array.from(grok.shell.tableViews).find((v) => !before.has(v.id))) !== undefined,
        'no table view opened');
      expect(opened!.dataFrame.rowCount > 0, true);
    } finally {
      opened?.close();
    }
  });

  test('an object the platform has no default action for is left alone', async () => {
    const views = Array.from(grok.shell.views).length;
    DG.ObjectHandler.runDefaultAction(42);
    expect(Array.from(grok.shell.views).length, views);
  });
});
