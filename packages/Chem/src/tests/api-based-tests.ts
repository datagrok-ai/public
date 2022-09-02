import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, test, delay, before, after} from '@datagrok-libraries/utils/src/test';
import {_testSearchSubstructure, _testSearchSubstructureAllParameters} from './utils';
import {_testFindSimilar, _testGetSimilarities} from './menu-tests-similarity-diversity';
import {testCsv, testSubstructure} from './substructure-search-tests';
import {getHTMLElementbyInnerText, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkHTMLElementbyInnerText, isColumnPresent} from './gui-utils';
import {_importSdf} from '../open-chem/sdf-importer';

category('server features', () => {
  test('testDescriptors', async () => {
    const t = grok.data.demo.molecules();
    await grok.chem.descriptors(t, 'smiles', ['MolWt', 'Lipinski']);
  });
});

category('chem exported', () => {
  test('findSimilar.api.sar-small', async () => {
    await _testFindSimilar(grok.chem.findSimilar);
  });

  test('getSimilarities.api.molecules', async () => {
    await _testGetSimilarities(grok.chem.getSimilarities);
  });

  test('substructureSearch', async () => {
    const df = DG.DataFrame.fromCsv(testCsv);
    const trueIndices = [0, 2];
    const bitset: DG.BitSet = await grok.chem.searchSubstructure(df.col('smiles')!, testSubstructure);
    const bitsetString = bitset.toBinaryString();
    const bitsetArray = [...bitsetString];
    for (let k = 0; k < trueIndices.length; k++) {
      expect(bitsetArray[trueIndices[k]] === '1', true);
      bitsetArray[trueIndices[k]] = '0';
    }
  });

  test('mcs', async () => {
    const t = DG.DataFrame.fromCsv(`smiles
      O=C1CN=C(c2ccccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3
      O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4
      O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3
      CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3`);
    const mcs = await grok.chem.mcs(t.col('smiles')!);
    expect(mcs, 'O=C1CN=C(C2CCCCC2)C2:C:C:C:C:C:2N1');
  });
});

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
      if (mapIdentifiersAction.innerText == 'Chem | Map Identifiers...')
          break;
      }
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

category('Chem: Similarity Search', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
  });

  test('chem.similaritySearch', async () => {
    
    smiles.currentRowIdx = 0;

    grok.shell.topMenu.find('Chem').find('Similarity Search...').click(); await delay(2000);

    isViewerPresent(Array.from(v.viewers), 'SimilaritySearchViewer');

    let structuresInViewer = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0].getElementsByClassName('d4-flex-wrap ui-div')[0] as HTMLElement;
    if (structuresInViewer.childElementCount != 10)
        throw 'molecules number inside the viewer is different than expected'; 

    let similarityViewer:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'SimilaritySearchViewer') {
            similarityViewer = Array.from(v.viewers)[i];
            break;
        }
    }

    similarityViewer!.props.distanceMetric = 'Dice';
    similarityViewer!.props.limit = 5    

    let similarityLable = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0].getElementsByClassName('ui-label')[0] as HTMLElement;
    if (similarityLable.innerText != '0.22')
        throw 'Expected Similarity Lable for 2nd molecule does not match the "Dice" metric';
    
        structuresInViewer = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0].getElementsByClassName('d4-flex-wrap ui-div')[0] as HTMLElement;
    if (structuresInViewer.childElementCount != 5)
        throw 'molecules number inside Similarity viewer is different than expected after change "Limit" property'; 

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'SimilaritySearchViewer') {
            throw 'SimilaritySearch viewer was not closed'
        }
    }
  });  
  after(async () => {
    grok.shell.closeAll();
  }); 
});