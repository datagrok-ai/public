import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';

category('Grid: Data Search', () => {
  let v: DG.TableView;
  let demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(await demog);
  });

  test('grid.dataSearch', async () => {
    let searchTab:HTMLElement | undefined;
    let tab;
    for (let i=0; i<document.getElementsByClassName('d4-accordion-pane-header').length; i++) {
        tab = document.getElementsByClassName('d4-accordion-pane-header')[i] as HTMLElement;
        if (tab.innerText == 'Search'){
            searchTab = tab;
            break;
        }
    }

    searchTab!.click();

    let searchInput:HTMLInputElement | undefined;
    let input;
    for (let i=0; i<document.getElementsByClassName('ui-input-editor').length; i++) {
        input = document.getElementsByClassName('ui-input-editor')[i] as HTMLInputElement;
        if (input.placeholder == 'Search'){
            searchInput = input;
            break;
        }
    }

    searchInput!.value = 'Asian'; await delay(500);

    let searchOptionsBtn = document.getElementsByClassName('d4-flex-row d4-flew-nowrap d4-search')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    searchOptionsBtn.click(); await delay(500);

    let filterMatchingAction:HTMLElement | undefined;
    let actionItem;
    for (let i=0; i<document.getElementsByClassName('d4-menu-item-label').length; i++) {
        actionItem = document.getElementsByClassName('d4-menu-item-label')[i] as HTMLElement;
        if (actionItem.innerText == 'Filter matching'){
            filterMatchingAction = actionItem;
            break;
        }
    }
    
    filterMatchingAction!.click(); await delay(500);

    if (demog.filter.trueCount != 762)
    throw 'Data searching error'    
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
