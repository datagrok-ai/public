import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Inputs', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.inputs', async () => {
    grok.shell.topMenu.find('Chem').find('Mutate...').click(); await delay(500);
    isDialogPresent('Mutate');

    expect(returnDialog('Mutate')!.input('Smiles').stringValue, 'CN1C(CC(O)C1=O)C1=CN=CC=C1');

    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); await delay(2000);

    expect(grok.shell.t.name, 'mutations')
  });
    
  after(async () => {
    grok.shell.closeAll();
  });  
});
