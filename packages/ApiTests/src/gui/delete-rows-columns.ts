import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';

category('Grid: Delete Rows Columns', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(await demog);
  });

  test('grid.deleteRows', async () => {
    grok.shell.topMenu.find('Select').find('Random...').click(); await delay(1000);
    isDialogPresent('Select Random Rows')
  
    let input25 = document.getElementsByName('label-25%')[0] as HTMLElement;
    input25!.click(); await delay(500);

    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
    okButton!.click(); await delay(1000);

    let removeRowsBtn = document.getElementsByClassName('svg-remove-selected-rows')[0] as HTMLElement;
    removeRowsBtn.click(); await delay(500);

    if (demog.rowCount == 1000)
        throw 'rows are not deleted';
  });

  test('grid.deleteCols', async () => {
    let demogHeader = document.getElementsByClassName('d4-ribbon-name')[0] as HTMLElement;
    demogHeader.click(); await delay(500);

    let columnsSectionOnPP = document.getElementsByName('div-section--Columns')[0] as HTMLElement;
    columnsSectionOnPP.click(); await delay(500);

    let diseaseColOnPP = document.getElementsByName('span-disease')[0] as HTMLElement;
    diseaseColOnPP.click(); await delay(500);

    let actionsSectionOnPP = document.getElementsByName('div-section--Actions')[0] as HTMLElement;
    actionsSectionOnPP.click(); await delay(500);

    //here we find Remove link action PP
    let removeLinkAction:HTMLElement | undefined;
    let linkActionLabel;
    for (let i=0; i<document.getElementsByClassName('d4-link-action').length; i++) {
        linkActionLabel = document.getElementsByClassName('d4-link-action')[i] as HTMLElement;
        if (linkActionLabel.innerText == 'Remove'){
            removeLinkAction = linkActionLabel;
            break;
        }
    }

    removeLinkAction!.click(); await delay(500);   
    
    if (demog.columns.byName('disease') != null)
      throw 'disease column was not deleted'
  }); 

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
