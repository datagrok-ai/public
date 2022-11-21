import {after, before, category, test, delay} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {waitForElement} from '../gui/gui-utils';


category('Connections', () => {
  let v: DG.ViewBase;
  let table: HTMLElement;

  before(async () => {
    const mng: DG.TabPane = grok.shell.sidebar.getPane('Manage');
    await (mng.content.querySelector('[data-view=connections]') as HTMLElement).click();
    const nw = await waitForElement('[data-link="/e/ApiTests.PostgreSQLNorthwind"]',
      'cannot find Northwind connection');
    v = grok.shell.v;
    nw.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await delay(50);
    const bs = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Browse schema')) as HTMLElement;
    bs.click();
    v.close();
    table = await waitForElement('.d4-sketch-item', 'cannot find table');
    v = grok.shell.v;
  });

  test('getAll', async () => {
    table.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await delay(50);
    const getAll = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Get All')) as HTMLElement;
    getAll.click();
    await waitForElement('canvas', 'Get All does not work', 5000);
    grok.shell.v.close();
  });

  test('getTop100', async () => {
    table.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true}));
    await delay(50);
    const getTop100 = Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((el) => (el as HTMLElement).innerText.includes('Get Top 100')) as HTMLElement;
    getTop100.click();
    await waitForElement('canvas', 'Get Top 100 does not work', 5000);
    grok.shell.v.close();
  });

  after(async () => {
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });
});
