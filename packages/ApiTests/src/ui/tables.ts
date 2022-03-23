import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Tables', () => {
  let v: DG.View;

  let arr = [
    { key: 'First', value: true },
    { key: 'Second', value: false },
    { key: 'Third', value: false },
  ];

  let table = ui.table(arr, (item, idx) => [item.key, item.value]);
  let tableFromMap = ui.tableFromMap(arr);
  let HTMLTable = DG.HtmlTable.create(arr, (item: any, idx: any) => [item.key, item.value]);

  before(async () => {
    v = grok.shell.newView('');
  });

  test('table.root', async () => {
    checkHTMLElement('Table', table, v, 'table');
  });

  test('tableFromMap.root', async () => {
    checkHTMLElement('tableFromMap', tableFromMap, v, 'table');
  });

  test('HTMLTable.root', async () => {
    checkHTMLElement('HTMLTable', HTMLTable.root, v, 'table');
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });

});