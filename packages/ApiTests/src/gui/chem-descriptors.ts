import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Descriptors', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.descriptors.ui', async () => {
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
    let mapIdentifiersAction:HTMLElement;
    for (let i = 0; i < actions.length; i++ ) {        
        mapIdentifiersAction = actions[i] as HTMLElement;
        if (mapIdentifiersAction.innerText == 'Chem | Descriptors...')
            break;
    }
    mapIdentifiersAction!.click();

    await delay(500);    
    isDialogPresent('Descriptors')

    let discriptorsGroups = returnDialog('Descriptors')!.root.getElementsByClassName('d4-tree-view-node');
    let lipinskiGroup:HTMLElement | undefined = undefined;
    for (let i = 0; i < discriptorsGroups.length; i++ ) {  
        lipinskiGroup = discriptorsGroups[i] as HTMLElement;
        if (lipinskiGroup.innerText == 'Lipinski')
            break;
    }

    let lipinskiChekbox:HTMLElement | undefined;
    for (let i = 0; i < lipinskiGroup!.childNodes.length; i++ ) {   
        lipinskiChekbox = lipinskiGroup!.childNodes[i] as HTMLElement;
        if (lipinskiChekbox.tagName == 'INPUT')
            break;
    }

    lipinskiChekbox!.click();

    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); await delay(500);

    isColumnPresent(smiles.columns, 'FractionCSP3');
    isColumnPresent(smiles.columns, 'NumAromaticCarbocycles');
    isColumnPresent(smiles.columns, 'NumHAcceptors');
    isColumnPresent(smiles.columns, 'NumHeteroatoms');
    isColumnPresent(smiles.columns, 'NumRotatableBonds');
    isColumnPresent(smiles.columns, 'RingCount');

    v.close();
  });
  test('chem.descriptors.api', async () => {

    grok.chem.descriptors(grok.data.testData('molecules', 100), 'smiles', ['MolWt', 'Lipinski'])
    .then(function (table) {
      grok.shell.addTableView(table);
    });

    await delay(1000);

    isColumnPresent(grok.shell.t.columns, 'MolWt');
    isColumnPresent(grok.shell.t.columns, 'NumAromaticCarbocycles');
    isColumnPresent(grok.shell.t.columns, 'NumHAcceptors');
    isColumnPresent(grok.shell.t.columns, 'NumHeteroatoms');
    isColumnPresent(grok.shell.t.columns, 'NumRotatableBonds');
    isColumnPresent(grok.shell.t.columns, 'RingCount');
  });  
    
   after(async () => {
    grok.shell.closeAll();
  });   
});
