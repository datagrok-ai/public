import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: R Groups Analysis', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
  });

  test('rGroupsAnalysis.ui', async () => {    
    grok.shell.topMenu.find('Chem').find('R-Groups Analysis......').root.click(); await delay(2000);

    isDialogPresent('R-Groups Analysis');

    let mcsBtn = document.getElementsByClassName('d4-dialog-contents dlg-r-groups-analysis ui-form ui-panel')[0].getElementsByTagName('button')[0] as HTMLElement;
    mcsBtn.click(); await delay(7000);

    let prefixInput = document.getElementsByClassName('d4-dialog-contents dlg-r-groups-analysis ui-form ui-panel')[0].getElementsByClassName('ui-input-text ui-input-root')[0].getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;
    prefixInput.value = 'UITestRGroup';

    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); await delay(5000);

    isViewerPresent(Array.from(v.viewers), 'Trellis plot');
    isColumnPresent(smiles.columns, 'UITestRGroup1');
    isColumnPresent(smiles.columns, 'UITestRGroup2');
    isColumnPresent(smiles.columns, 'UITestRGroup3');
    isColumnPresent(smiles.columns, 'UITestRGroup4');

    }); 
  test('rGroupsAnalysis.api', async () => {
    v.resetLayout(); await delay(500);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Trellis plot')
          throw 'Table view was not cleared';
      }

    let smilesCol = smiles.col('smiles') as DG.Column;
    let mcs = await grok.chem.mcs(smilesCol);      
    await grok.chem.rGroup(smiles, 'smiles', mcs);
    
    isViewerPresent(Array.from(v.viewers), 'Trellis plot');
    isColumnPresent(smiles.columns, 'R1');
    isColumnPresent(smiles.columns, 'R2');
    isColumnPresent(smiles.columns, 'R3');
    isColumnPresent(smiles.columns, 'R4');
    
  });  

 
  after(async () => {
    grok.shell.closeAll();
  }); 
});
