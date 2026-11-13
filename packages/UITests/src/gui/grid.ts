import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

import { category, test, waitUntil, before } from '@datagrok-libraries/test/src/test';
import {getHTMLElementbyInnerText, showToolbox} from './gui-utils';
import { checkDialog } from './gui-utils';

category('GUI: Grid', () => {
  before(async () => {
    showToolbox();
    grok.shell.windows.showContextPanel = true;
  });

  test('grid.dataSearch', async () => {
    const demog = grok.data.demo.demog(100);
    const v = grok.shell.addTableView(demog);
    await waitUntil(() => grok.shell.v == v, 'table view is not current');
    const searchTab = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Search');
    if (!searchTab?.classList.contains('expanded'))
      searchTab!.click();

    let searchInput: HTMLInputElement | undefined;
    await waitUntil(() => {
      searchInput = Array.from(searchTab!.parentElement!.getElementsByClassName('ui-input-editor'))
        .find((el) => (el as HTMLInputElement).placeholder == 'Search') as HTMLInputElement | undefined;
      return searchInput != undefined;
    }, 'cannot find search input');
    searchInput!.value = 'other';
    searchInput!.dispatchEvent(new Event('input'));

    const searchOptionsBtn = document.getElementsByClassName('d4-flex-row d4-flew-nowrap d4-search')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    searchOptionsBtn.click();
    await waitUntil(() => getHTMLElementbyInnerText('d4-menu-item-label', 'Filter matching') != undefined,
      'cannot find Filter matching option');

    const filterMatchingAction = getHTMLElementbyInnerText('d4-menu-item-label', 'Filter matching');
    filterMatchingAction!.click();
    await waitUntil(() => demog.filter.trueCount == 8, 'error in Filter matching');
  });

  test('grid.deleteRows', async () => {
    const demog = grok.data.demo.demog(100);
    const v = grok.shell.addTableView(demog);
    await waitUntil(() => grok.shell.v == v, 'table view is not current');
    grok.shell.topMenu.find('Select').find('Random Rows...').click();
    await waitUntil(() => checkDialog('Select Random Rows'), 'cannot find dialog');

    const input25 = Array.from(document.querySelectorAll('.d4-link-label'))
      .find((el) => el.textContent === '25%') as HTMLElement;
    input25!.click();
    await waitUntil(() => demog.selection.trueCount == 25, 'random rows are not selected');
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await waitUntil(() => !checkDialog('Select Random Rows'), 'dialog is not closed');

    const removeRowsBtn = document.getElementsByClassName('svg-remove-selected-rows')[0] as HTMLElement;
    removeRowsBtn.click();
    await waitUntil(() => demog.rowCount == 75, 'rows are not deleted');
  });

  test('grid.deleteCols', async () => {
    const demog = grok.data.demo.demog(10);
    demog.currentCell = demog.cell(0, 'disease');
    grok.shell.addTableView(demog);
    await waitUntil(() => document.querySelector('canvas') !== null, 'cannot load table');
    grok.shell.o = demog.getCol('disease');

    await waitUntil(() => Array.from(document.querySelectorAll('.d4-accordion-title'))
      .find((el) => el.textContent == 'disease') !== undefined, 'cannot load column context panel');

    const actionsSectionOnPP = Array.from(document.getElementsByClassName('grok-entity-prop-panel')[0]
      .querySelectorAll('.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Actions') as HTMLElement;
    if (Array.from(actionsSectionOnPP.classList).find((el) => el == 'expanded') == undefined)
      actionsSectionOnPP.click();

    await waitUntil(() => actionsSectionOnPP.classList.contains('expanded'), 'Actions pane is not expanded');

    const removeLinkAction = Array.from(document.querySelectorAll('.d4-link-action'))
      .find((el) => el.textContent === 'Remove') as HTMLElement;
    removeLinkAction!.click();
    await waitUntil(() => demog.columns.byName('disease') == null, 'disease column was not deleted');
  });

  test('grid.filters', async () => {
    const demog = grok.data.demo.demog(100);
    grok.shell.addTableView(demog);
    await waitUntil(() => document.querySelector('canvas') !== null, 'cannot load table');

    await waitUntil(() => grok.shell.topMenu.find('Select') != null, 'cannot find Select in top menu');
    grok.shell.topMenu.find('Select').find('Random Rows...').click();
    await waitUntil(() => checkDialog('Select Random Rows'), 'cannot find dialog');

    await waitUntil(() => Array.from(document.querySelectorAll('.d4-link-label'))
      .find((el) => el.textContent === '25%') !== undefined, 'cannot find 25% in the dialog');

    const input25 = Array.from(document.querySelectorAll('.d4-link-label'))
      .find((el) => el.textContent === '25%') as HTMLElement;
    input25!.click();
    await waitUntil(() => demog.selection.trueCount == 25, 'random rows are not selected');

    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await waitUntil(() => !checkDialog('Select Random Rows'), 'dialog is not closed');

    const cp = document.querySelector('.grok-prop-panel') as HTMLElement;
    await waitUntil(() => Array.from(cp.querySelectorAll('.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Actions') !== undefined, 'cannot find Actions on context panel');
    const actionsSectionOnPP = Array.from(cp.querySelectorAll('.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Actions') as HTMLElement;
    if (!actionsSectionOnPP?.classList.contains('expanded'))
      actionsSectionOnPP.click();
    let filterLinkAction: HTMLElement | undefined;
    await waitUntil(() => {
      filterLinkAction = Array.from(document.querySelectorAll('.d4-link-action'))
        .find((el) => el.textContent === 'Filter Rows') as HTMLElement | undefined;
      return filterLinkAction !== undefined;
    }, 'cannot find Filter Rows action');

    filterLinkAction!.click();
    await waitUntil(() => demog.filter.trueCount == 25, 'Error in filtering');
  });
}, {owner: 'dkovalyov@datagrok.ai'});
