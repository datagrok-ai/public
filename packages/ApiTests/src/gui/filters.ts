import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';

category('Grid: Filters', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(await demog);
  });

  test('grid.filters', async () => {
    grok.shell.topMenu.find('Select').find('Random...').root.click(); await delay(1000);
    isDialogPresent('Select Random Rows')
  
    let input25 = document.getElementsByName('label-25%')[0] as HTMLElement;
    input25!.click(); await delay(500);

    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
    okButton!.click(); await delay(1000);

    let actionsSectionOnPP = document.getElementsByName('div-section--Actions')[0] as HTMLElement;
    actionsSectionOnPP.click(); await delay(500);

    let filterLinkAction:HTMLElement | undefined;
    let linkActionLabel;
    for (let i=0; i<document.getElementsByClassName('d4-link-action').length; i++) {
        linkActionLabel = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
        if (linkActionLabel.innerText == 'Filter Rows')
        filterLinkAction = linkActionLabel;
        break;
    }

    filterLinkAction!.click(); await delay(1000);

    if (demog.filter.trueCount != 250)
     throw 'Error in filtering'
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
