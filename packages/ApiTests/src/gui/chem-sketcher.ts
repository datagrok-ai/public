import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Sketcher', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.sketcher', async () => {
  /*   grok.shell.topMenu.find('Chem').find('Sketcher').click(); await delay(500);

    isDialogPresent('Sketcher');

    returnDialog('Sketcher')!.close(); 
    await delay(3000); */
    smiles.currentRowIdx = 5; await delay(1000);

    let panels = document.getElementsByClassName('grok-prop-panel')[0].getElementsByClassName('d4-accordion-pane-header');
    let actionsPanel:HTMLElement;
    for (let i = 0; i < panels.length; i++ ){        
      actionsPanel = panels[i] as HTMLElement;
      if (actionsPanel.innerText == 'Actions')
          break;}
    actionsPanel!.click(); await delay(500);

    let sketchBtn = getHTMLElementbyInnerText('d4-link-action', 'Sketch')
    sketchBtn!.click(); await delay(1000);

    isDialogPresent('Sketcher');

  }); 

  after(async () => {
    v.close();
    grok.shell.closeAll();
  }); 
});
