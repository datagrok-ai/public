import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {getHTMLElementbyInnerText} from './gui-utils';
import {checkDialog} from './gui-utils';

category('GUI: Grid', () => {
  test('grid.dataSearch', async () => {
    const demog = grok.data.demo.demog(1000);
    const v = grok.shell.addTableView(demog);
    await awaitCheck(() => grok.shell.v == v);

    try {
      const searchTab = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Search');
      searchTab!.click();

      let searchInput:HTMLInputElement | undefined;
      let input;
      for (let i=0; i<document.getElementsByClassName('ui-input-editor').length; i++) {
        input = document.getElementsByClassName('ui-input-editor')[i] as HTMLInputElement;
        if (input.placeholder == 'Search') {
          searchInput = input;
          break;
        }
      }

      await awaitCheck(() => searchInput != undefined);
      searchInput!.value = 'Asian'; await delay(100);
      searchInput!.dispatchEvent(new Event('input'));

      const searchOptionsBtn = document.getElementsByClassName('d4-flex-row d4-flew-nowrap d4-search')[0]
        .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
      searchOptionsBtn.click();

      await awaitCheck(() => getHTMLElementbyInnerText('d4-menu-item-label', 'Filter matching') != undefined);

      const filterMatchingAction = getHTMLElementbyInnerText('d4-menu-item-label', 'Filter matching');
      filterMatchingAction!.click();

      await awaitCheck(() => demog.filter.trueCount == 762);
    } finally {
      grok.shell.closeAll();
    }
  });

  test('grid.deleteRows', async () => {
    const demog = grok.data.demo.demog(1000);
    const v = grok.shell.addTableView(demog);
    await awaitCheck(() => grok.shell.v == v);

    try {
      grok.shell.topMenu.find('Select').find('Random...').click();
      await awaitCheck(() => checkDialog('Select Random Rows'));

      const input25 = Array.from(document.querySelectorAll('.d4-link-label'))
        .find((el) => el.textContent === '25%') as HTMLElement;
      input25!.click();
      await awaitCheck(() => demog.selection.trueCount == 250);

      const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
      okButton.click();
      await delay(200);

      const removeRowsBtn = document.getElementsByClassName('svg-remove-selected-rows')[0] as HTMLElement;
      removeRowsBtn.click();

      await awaitCheck(() => demog.rowCount == 750);

      if (demog.rowCount == 1000)
        throw new Error('rows are not deleted');
    } finally {
      grok.shell.closeAll();
    }
  });

  test('grid.deleteCols', async () => {
    const demog = grok.data.demo.demog(10);
    demog.currentCell = demog.cell(0, 'disease');
    const v = grok.shell.addTableView(demog);
    await awaitCheck(() => grok.shell.v == v);
    await delay(1000);
    grok.shell.o = demog.getCol('disease');

    try {
      await awaitCheck(() => Array.from(document.querySelectorAll('.d4-accordion-title'))
        .find((el) => el.textContent == 'disease') !== undefined, 'cannot load column context panel', 3000);

      const actionsSectionOnPP = Array.from(document.getElementsByClassName('grok-entity-prop-panel')[0]
        .querySelectorAll('.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Actions') as HTMLElement;
      if (Array.from(actionsSectionOnPP.classList).find((el) => el == 'expanded') == undefined)
        actionsSectionOnPP.click();

      await awaitCheck(() => Array.from(actionsSectionOnPP.classList)
        .find((el) => el == 'expanded') !== undefined);

      const removeLinkAction = Array.from(document.querySelectorAll('.d4-link-action'))
        .find((el) => el.textContent === 'Remove') as HTMLElement;
      removeLinkAction!.click();

      await awaitCheck(() => demog.columns.byName('disease') == null, 'column is not deleted', 2000);

      if (demog.columns.byName('disease') != null)
        throw new Error('disease column was not deleted');
    } finally {
      grok.shell.closeAll();
    }
  });

  test('grid.filters', async () => {
    const demog = grok.data.demo.demog(1000);
    const v = grok.shell.addTableView(demog);
    await awaitCheck(() => grok.shell.v == v);

    try {
      grok.shell.topMenu.find('Select').find('Random...').click();
      await awaitCheck(() => checkDialog('Select Random Rows'));

      const input25 = Array.from(document.querySelectorAll('.d4-link-label'))
        .find((el) => el.textContent === '25%') as HTMLElement;
      input25!.click();
      await awaitCheck(() => demog.selection.trueCount == 250);

      const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
      okButton.click();
      await delay(200);

      const actionsSectionOnPP = Array.from(document.getElementsByClassName('grok-prop-panel')[0]
        .querySelectorAll('.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Actions') as HTMLElement;
      if (Array.from(actionsSectionOnPP.classList).find((e) => e == 'expanded') == undefined)
        actionsSectionOnPP.click();

      await awaitCheck(() => Array.from(actionsSectionOnPP.classList).find((e) => e == 'expanded') !== undefined);

      const filterLinkAction = Array.from(document.querySelectorAll('.d4-link-action'))
        .find((el) => el.textContent === 'Filter Rows') as HTMLElement;

      filterLinkAction!.click();
      await awaitCheck(() => demog.filter.trueCount == 250);

      if (demog.filter.trueCount != 250)
        throw new Error('Error in filtering');
    } finally {
      grok.shell.closeAll();
    }
  });

  after(async () => {
    grok.shell.closeAll();
  });
});
