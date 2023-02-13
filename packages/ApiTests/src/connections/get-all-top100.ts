import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Connections', () => {
  let schema: DG.ViewBase;
  let table: HTMLElement;
  let nw: HTMLElement;

  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Manage');
    await (mng.content.querySelector('[data-view=connections]') as HTMLElement).click();
    await awaitCheck(() => document.querySelector('[data-link="/e/ApiTests.PostgresTest"]') !== null,
      'cannot find Northwind connection', 3000);
    nw = document.querySelector('[data-link="/e/ApiTests.PostgresTest"]') as HTMLElement;
  });

  test('getAll', async () => {
    await browseSchema();
    table.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await awaitCheck(() => document.querySelector('.d4-menu-popup') !== null, 'cannot find context menu', 2000);
    const getAll = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Get All')) as HTMLElement;
    getAll.click();
    await awaitCheck(() => document.querySelector('canvas') !== null, 'Get All does not work', 5000);
    grok.shell.v.close();
    schema.close();
  });

  test('getTop100', async () => {
    await browseSchema();
    table.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await awaitCheck(() => document.querySelector('.d4-menu-popup') !== null, 'cannot find context menu', 2000);
    const getTop100 = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Get Top 100')) as HTMLElement;
    getTop100.click();
    await awaitCheck(() => document.querySelector('canvas') !== null, 'Get Top 100 does not work', 5000);
    grok.shell.v.close();
    schema.close();
  });

  after(async () => {
    grok.shell.closeAll();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });

  async function browseSchema() {
    nw.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await awaitCheck(() => document.querySelector('.d4-menu-popup') !== null, 'cannot find context menu');
    const bs = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Browse schema')) as HTMLElement;
    bs.click();
    await awaitCheck(() => document.querySelector('.d4-sketch-item') !== null, 'cannot find table', 10000);
    table = document.querySelector('.d4-sketch-item') as HTMLElement;
    schema = grok.shell.v;
  }
});
