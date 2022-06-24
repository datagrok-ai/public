import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';

category('Dialog: SPE', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(await demog);
  });

  test('dialogs.spe', async () => {
    grok.shell.topMenu.find('Tools').find('Data Science').find('Stochastic Proximity Embedding...').root.click(); await delay(1000);

    const okButton = document.getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
        okButton!.click(); await delay(2000);

    //WIP, SPE Feature currently not working (GROK-10250)


  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
