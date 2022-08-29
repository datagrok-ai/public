import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Map Identifiers', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.mapIdentifiers', async () => {
    let smilesCol = smiles.columns.byName("smiles");
    grok.shell.o = smilesCol;

    await delay(1000);

    let panels = document.getElementsByClassName('grok-prop-panel')[0].getElementsByClassName('d4-accordion-pane-header');
    let actionsPanel:HTMLElement;
    for (let i = 0; i < panels.length; i++ ){        
      actionsPanel = panels[i] as HTMLElement;
      console.log('PANEL: ' + actionsPanel.innerText);
      if (actionsPanel.innerText == 'Actions')
          break;}
    actionsPanel!.click(); 
    
    await delay(500);

    let actions = document.getElementsByClassName('grok-prop-panel')[0].getElementsByClassName('d4-link-action');
    let mapIdentifiersAction:HTMLElement;
    for (let i = 0; i < actions.length; i++ ){        
      mapIdentifiersAction = actions[i] as HTMLElement;
      if (mapIdentifiersAction.innerText == 'Chem | Map Identifiers...')
          break;}
    mapIdentifiersAction!.click();

    await delay(500);    
    isDialogPresent('Chem Map Identifiers')

    setDialogInputValue('Chem Map Identifiers', 'To Source', 'mcule');
    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); await delay(3000);
    
    isColumnPresent(smiles.columns, 'mcule');
  });
    
  after(async () => {
    v.close();
    grok.shell.closeAll();
  });  
});
