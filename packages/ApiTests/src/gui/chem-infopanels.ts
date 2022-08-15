import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Info Panels', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
    await delay(3000);
  });

  test('chem.infoPanel.gasteigerPartialCharges', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Gasteiger Partial Charges');

    let gastaigerPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Gasteiger Partial Charges')
    gastaigerPanel!.click(); await delay(4000);

    if (document.getElementsByClassName('grok-scripting-image-container-info-panel').length != 1)
        throw 'script output was not rendered in the panel'

    let pecilIcon = document.getElementsByClassName('grok-icon fal fa-pencil')[0] as HTMLElement;
    pecilIcon!.click(); await delay(200);

    let contours = document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-gasteiger_partial_charges')[0].getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;
    contours.value = '15';

    let applyBtn = document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-gasteiger_partial_charges')[0].getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
    applyBtn!.click(); await delay(2000);

    gastaigerPanel!.click(); await delay(500);
  }); 
  test('chem.infoPanel.identifiers', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Identifiers');

    let identifiersPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Identifiers')
    identifiersPanel!.click(); await delay(2000);

    checkHTMLElementbyInnerText('ui-link d4-link-external', 'DTXSID00639834');
    checkHTMLElementbyInnerText('ui-link d4-link-external', '24204979');
    checkHTMLElementbyInnerText('ui-link d4-link-external', 'CHEMBL2262434');
        
    identifiersPanel!.click(); await delay(500);
  }); 
  test('chem.infoPanel.structure2D', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Structure 2D');

    let structure2DPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Structure 2D')
    structure2DPanel!.click(); await delay(1000);

    if (document.getElementsByClassName('d4-flex-col ui-div chem-mol-box').length != 1)
        throw 'canvas with structure was not rendered in the panel'    
        
    structure2DPanel!.click(); await delay(500);
  }); 
  test('chem.infoPanel.structure3D', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Structure 3D');

    let structure3DPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Structure 3D')
    structure3DPanel!.click(); await delay(2000);

    if (document.getElementsByClassName('d4-ngl-viewer ui-div').length != 1)
        throw 'canvas with 3D structure was not rendered in the panel'    
        
    structure3DPanel!.click(); await delay(500);
  }); 
  test('chem.infoPanel.properties', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Properties');

    let propertiesPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Properties')
    propertiesPanel!.click(); await delay(2000);

    if (document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-properties')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'table with properties was not rendered in the panel'    
        
    propertiesPanel!.click(); await delay(500);
  }); 
  test('chem.infoPanel.toxicity', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Toxicity');

    let toxicityPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Toxicity')
    toxicityPanel!.click(); await delay(2000);

    if (document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-toxicity')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'table with toxicity was not rendered in the panel'    
        
    toxicityPanel!.click(); await delay(500);
  });  
  test('chem.infoPanel.drugLikeness', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Drug Likeness');

    let drugLikenessPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Drug Likeness')
    drugLikenessPanel!.click(); await delay(2000);

    if (document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-drug_likeness')[0].getElementsByClassName('d4-flex-col ui-div').length != 19)
        throw 'number of displayed canvases with molecules does not match the expected'    
        
        drugLikenessPanel!.click(); await delay(500);
  });  
  test('chem.infoPanel.structuralAlerts', async () => {
    smiles.currentRowIdx = 5; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Structural Alerts');

    let structuralAlertsPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Structural Alerts')
    structuralAlertsPanel!.click(); await delay(2000);

    if (document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-structural_alerts')[0].getElementsByClassName('d4-flex-col ui-div').length != 10)
        throw 'number of displayed canvases with alerts does not match the expected'    
        
        structuralAlertsPanel!.click(); await delay(500);
  });  
  after(async () => {
    v.close();
    grok.shell.closeAll();
  }); 
});
