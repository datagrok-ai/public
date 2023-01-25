import {after, before, category, test, delay, awaitCheck} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('Connections', () => {
  let v: DG.ViewBase;
  let table: HTMLElement;

  before(async () => {
    try {
      const mng: DG.TabPane = grok.shell.sidebar.getPane('Manage');
      await (mng.content.querySelector('[data-view=connections]') as HTMLElement).click();
      await awaitCheck(() => document.querySelector('[data-link="/e/ApiTests.PostgreSQLNorthwind"]') !== null,
        'cannot find Northwind connection', 3000);
      const nw = document.querySelector('[data-link="/e/ApiTests.PostgreSQLNorthwind"]');
      v = grok.shell.v;
      nw!.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
      await delay(50);
      const bs = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((el) => (el as HTMLElement).innerText.includes('Browse schema')) as HTMLElement;
      bs.click();
      v.close();
      await awaitCheck(() => document.querySelector('.d4-sketch-item') !== null, 'cannot find table', 3000);
      table = document.querySelector('.d4-sketch-item') as HTMLElement;
      v = grok.shell.v;
    } catch (e) {
      console.log(e);
    }
  });

  test('getAll', async () => {
    table.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await delay(50);
    const getAll = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Get All')) as HTMLElement;
    getAll.click();
    await awaitCheck(() => document.querySelector('canvas') !== null, 'Get All does not work', 5000);
    grok.shell.v.close();
  });

  test('getTop100', async () => {
    table.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await delay(50);
    const getTop100 = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Get Top 100')) as HTMLElement;
    getTop100.click();
    await awaitCheck(() => document.querySelector('canvas') !== null, 'Get Top 100 does not work', 5000);
    grok.shell.v.close();
  });

  after(async () => {
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });
});
