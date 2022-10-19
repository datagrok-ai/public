import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Fingerprints', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.fingerprints', async () => {
    let smilesCol = smiles.columns.byName("smiles");
    grok.shell.o = smilesCol;

    await delay(1000);

    let panels = document.getElementsByClassName('grok-prop-panel')[0].getElementsByClassName('d4-accordion-pane-header');
    let actionsPanel:HTMLElement;
    for (let i = 0; i < panels.length; i++ ) {        
      actionsPanel = panels[i] as HTMLElement;
      if (actionsPanel.innerText == 'Actions')
          break;
      }
    actionsPanel!.click(); 
    
    await delay(500);

    let actions = document.getElementsByClassName('grok-prop-panel')[0].getElementsByClassName('d4-link-action');
    let fingerprintsAction:HTMLElement;
    for (let i = 0; i < actions.length; i++ ) {        
        fingerprintsAction = actions[i] as HTMLElement;
      if (fingerprintsAction.innerText == 'Chem | Fingerprints...')
          break;
      }
      fingerprintsAction!.click();

    await delay(500);    
    isDialogPresent('Fingerprints')

    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); await delay(3000);
    
    isColumnPresent(smiles.columns, 'Fingerprints');
  });
    
  after(async () => {
    v.close();
    grok.shell.closeAll();
  });  
});
