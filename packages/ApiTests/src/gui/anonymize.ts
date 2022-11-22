import {after, before, awaitCheck ,category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, waitForElement} from './gui-utils';

category('GUI: Anonymize', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('dialogs.anonymize', async () => {
    grok.shell.topMenu.find('Data').find('Anonymize...').click(); 

    function checkDialog(dialogTitle:string):boolean {
      let check = false;
      for (let i=0; i < DG.Dialog.getOpenDialogs().length; i++) {
        if (DG.Dialog.getOpenDialogs()[i].title == dialogTitle) {
          check = true;
          break;
        }
      }
      return check;  
    }

    await awaitCheck(() => {return checkDialog('Anonymize Data')});

    setDialogInputValue('Anonymize Data', 'Number randomization factor', 1);

    let okButton:HTMLElement | undefined;
    let button;
    for (let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++) {
      button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
      if (button.innerText == 'OK')
        okButton = button;
    }
    okButton!.click();
    await awaitCheck(() => {return !checkDialog('Anonimize Data')});
  });

  after(async () => {
    grok.shell.closeAll();
  });
});
