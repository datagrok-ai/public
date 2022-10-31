import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';

category('Chem: ToInchi', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.toInchi', async () => {
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
    let toInchiAction:HTMLElement;
    for (let i = 0; i < actions.length; i++ ) {        
        toInchiAction = actions[i] as HTMLElement;
      if (toInchiAction.innerText == 'Chem | To InchI')
          break;
      }
    toInchiAction!.click();

    await delay(1500);    
        
    isColumnPresent(smiles.columns, 'inchi');
  });
    
  after(async () => {
    v.close();
    grok.shell.closeAll();
  });  
});
