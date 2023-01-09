import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, after, expect, test, delay, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {setDialogInputValue, isColumnPresent} from './gui-utils';


category('UI', () => {
  let v: DG.TableView;
  let smiles: DG.DataFrame;

  test('similarity search', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.topMenu.find('Chem').find('Similarity Search...').click();
    await awaitCheck(() => document.querySelector('.d4-similaritysearchviewer') !== null, 'cannot load Similarity Search viewer', 2000);
    const similarityViewer = Array.from(v.viewers)[1];
    await awaitCheck(() => similarityViewer.root.querySelectorAll('.chem-canvas').length === 10,
      'molecules number inside Similarity viewer is different than expected', 3000);

    similarityViewer!.props.distanceMetric = 'Dice';
    similarityViewer!.props.limit = 5;    
    await awaitCheck(() => similarityViewer.root.querySelectorAll('.chem-canvas').length === 5,
      'molecules number inside Similarity viewer is different than expected after change "Limit" property', 3000);

    const similarityLable = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0]
      .getElementsByClassName('ui-label')[0] as HTMLElement;
    if (similarityLable.innerText != '0.22')
      throw 'Expected Similarity Lable for 2nd molecule does not match the "Dice" metric';
    
    // const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
    //   .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    // protpertiesBtn.click();
    // await awaitCheck(() => document.getElementsByClassName('property-grid-base property-grid-disable-selection').length > 0,
    //   'properties table does not open', 3000);

    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'SimilaritySearch viewer was not closed', 1000);
    v.close();
  });
  
  test('diversity search', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.topMenu.find('Chem').find('Diversity Search...').click();
    await awaitCheck(() => document.querySelector('.d4-diversitysearchviewer') !== null, 'cannot load Diversity Search viewer', 2000);
    const dsvRoot = document.querySelector('.d4-diversitysearchviewer') as HTMLElement;
    await awaitCheck(() => dsvRoot.querySelectorAll('.chem-canvas').length === 10, 'molecules number != 10', 3000);
    const dsv = Array.from(v.viewers)[1];
    dsv.setOptions({
      distanceMetric: 'Dice',
      size: 'normal',
    });
    dsv!.props.limit = 5;    
    await awaitCheck(() => dsvRoot.querySelectorAll('.chem-canvas').length === 5, 'molecules number != 5', 3000);
    v.close();
  });
  
  test('descriptors', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Toxicity') !== undefined;
    }, 'cannot find Toxicity property', 5000);
    const smilesCol = smiles.columns.byName('smiles');
    grok.shell.o = smilesCol;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Details') !== undefined;
    }, 'cannot load Smiles column properties', 5000);
    
    const actions = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Actions') as HTMLElement;
    if (!actions.classList.contains('expanded')) await actions.click();

    const descriptors = Array.from(pp.querySelectorAll('.d4-link-action'))
      .find((el) => el.textContent === 'Chem | Descriptors...') as HTMLElement;
    descriptors.click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length > 0, 'cannot find Descriptors dialog', 1000);
    const dialog = DG.Dialog.getOpenDialogs()[0].root;

    const lipinski = Array.from(dialog.querySelectorAll('.d4-tree-view-group-label'))
      .find((el) => el.textContent === 'Lipinski')!.previousSibling as HTMLElement;
    lipinski.click();

    const okButton = dialog.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => smiles.columns.length === 20, 'columns length != 20', 3000);

    isColumnPresent(smiles.columns, 'FractionCSP3');
    isColumnPresent(smiles.columns, 'NumAromaticCarbocycles');
    isColumnPresent(smiles.columns, 'NumHAcceptors');
    isColumnPresent(smiles.columns, 'NumHeteroatoms');
    isColumnPresent(smiles.columns, 'NumRotatableBonds');
    isColumnPresent(smiles.columns, 'RingCount');
    v.close();
  }, {skipReason: '#1183'});

  test('info panel: gasteiger', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
  
    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Gasteiger Partial Charges') !== undefined;
    }, 'cannot find Gasteiger Partial Charges property', 5000);
    
    const gpc = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Gasteiger Partial Charges') as HTMLElement;
    if (!gpc.classList.contains('expanded')) gpc.click();
    await awaitCheck(() => pp.querySelector('.grok-scripting-image-container-info-panel') !== null,
      'Gasteiger charges script output was not rendered in the panel', 10000);

    const pecilIcon = document.getElementsByClassName('grok-icon fal fa-pencil')[0] as HTMLElement;
    await pecilIcon!.click();
    const contours = document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-gasteiger_partial_charges')[0]
      .getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;
    contours.value = '15';
    const applyBtn = document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-gasteiger_partial_charges')[0]
      .getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
    applyBtn!.click();
    await delay(50);
    gpc.click();
    v.close();
  }, {skipReason: '#1183'});

  test('info panel: identifiers', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Identifiers') !== undefined;
    }, 'cannot find Identifiers property', 3000);

    const ih = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Identifiers') as HTMLElement;
    if (!ih.classList.contains('expanded')) ih.click();

    await awaitCheck(() => (ih.nextSibling as HTMLElement).querySelector('table') !== null, 'cannot load Identifiers', 3000);
    const it = ih.nextSibling as HTMLElement;
    for (const i of ['SCHEMBL5536145', '18722989', 'CHEMBL2262190']) {
      expect(Array.from(it.querySelectorAll('.ui-link.d4-link-external'))
        .find((el) => el.textContent === i) !== undefined, true);
    }
    ih.click(); await delay(10);
    v.close();
  }); 

  test('info panel: structure2D', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Structure 2D') !== undefined;
    }, 'cannot find Structure 2D property', 3000);

    const s2d = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Structure 2D') as HTMLElement;
    if (!s2d.classList.contains('expanded')) s2d.click();

    await awaitCheck(() => (s2d.nextSibling as HTMLElement).querySelector('.chem-canvas') !== null,
      'canvas with structure was not rendered in the panel', 3000);
    s2d.click(); await delay(10);
    v.close();
  });

  test('info panel: structure3D', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Structure 3D') !== undefined;
    }, 'cannot find Structure 3D property', 3000);

    const s3d = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Structure 3D') as HTMLElement;
    if (!s3d.classList.contains('expanded')) s3d.click();

    await awaitCheck(() => (s3d.nextSibling as HTMLElement).querySelector('canvas') !== null,
      'canvas with structure was not rendered in the panel', 10000);
    s3d.click(); await delay(10);
    v.close();
  });

  test('info panel: properties', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Properties') !== undefined;
    }, 'cannot find Properties property', 3000);

    const p = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Properties') as HTMLElement;
    if (!p.classList.contains('expanded')) p.click();

    await awaitCheck(() => (p.nextSibling as HTMLElement).querySelector('table') !== null,
      'table with properties was not rendered in the panel', 3000);
    p.click(); await delay(10);
    v.close();
  });

  test('info panel: toxicity', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Toxicity') !== undefined;
    }, 'cannot find Toxicity property', 3000);

    const t = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Toxicity') as HTMLElement;
    if (!t.classList.contains('expanded')) t.click();

    await awaitCheck(() => (t.nextSibling as HTMLElement).querySelector('table') !== null,
      'table with toxicity was not rendered in the panel', 3000);
    t.click(); await delay(10);
    v.close();
  });

  test('info panel: drug likeness', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Drug Likeness') !== undefined;
    }, 'cannot find Drug Likeness property', 3000);

    const dl = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Drug Likeness') as HTMLElement;
    if (!dl.classList.contains('expanded')) dl.click();

    await awaitCheck(() => (dl.nextSibling as HTMLElement).querySelectorAll('.d4-flex-col.ui-div').length === 50,
      'number of displayed canvases with molecules does not match the expected', 5000);
    dl.click(); await delay(10);
    v.close();
  });

  test('info panel: structural alerts', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    smiles.currentCell = smiles.cell(2, 'smiles');
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Structural Alerts') !== undefined;
    }, 'cannot find Structural Alerts property', 3000);

    const sa = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Structural Alerts') as HTMLElement;
    if (!sa.classList.contains('expanded')) sa.click();

    await awaitCheck(() => (sa.nextSibling as HTMLElement).querySelectorAll('.d4-flex-col.ui-div').length === 10,
      'number of displayed canvases with molecules does not match the expected', 5000);
    sa.click(); await delay(10);
    v.close();
  }, {skipReason: 'need to clear property panel'});

  test('chem inputs', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);

    grok.shell.topMenu.find('Chem').find('Mutate...').click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length > 0, 'cannot find Mutate dialog', 1000);
    const dialog = DG.Dialog.getOpenDialogs()[0];
    expect(dialog.input('Smiles').stringValue, 'CN1C(CC(O)C1=O)C1=CN=CC=C1');

    const okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.name === 'mutations', 'cannot find mutations table', 10000);
    await delay(10);
    grok.shell.v.close();
    grok.shell.closeTable(grok.shell.t);
    v.close();
  }, {skipReason: '#1183'});

  test('map identifiers', async () => {
    smiles = grok.data.demo.molecules(20);
    v = grok.shell.addTableView(smiles);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    const pp = document.querySelector('.grok-prop-panel') as HTMLElement;
    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Toxicity') !== undefined;
    }, 'cannot find Toxicity property', 5000);
    const smilesCol = smiles.columns.byName('smiles');
    grok.shell.o = smilesCol;

    await awaitCheck(() => {
      return Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
        .find((el) => el.textContent === 'Details') !== undefined;
    }, 'cannot load Smiles column properties', 5000);

    const actions = Array.from(pp.querySelectorAll('div.d4-accordion-pane-header'))
      .find((el) => el.textContent === 'Actions') as HTMLElement;
    if (!actions.classList.contains('expanded')) await actions.click();

    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'inchi');
    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('inchi'), 'cannot find inchi column', 10000);
 
    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'mcule');
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('mcule'), 'cannot find mcule column', 10000);

    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'chembl');
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('chembl'), 'cannot find chembl column', 10000);

    await callDialog();
    setDialogInputValue('Chem Map Identifiers', 'To Source', 'pubchem');
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click();
    await awaitCheck(() => grok.shell.t.columns.contains('pubchem'), 'cannot find pubchem column', 10000);
    v.close();

    async function callDialog() {
      const mi = Array.from(pp.querySelectorAll('.d4-link-action'))
        .find((el) => el.textContent === 'Chem | Map Identifiers...') as HTMLElement;
      mi.click();
      await awaitCheck(() => DG.Dialog.getOpenDialogs().length > 0, 'cannot find Chem Map Identifiers dialog', 1000);
    }
  }, {skipReason: 'need to clear property panel'});

  after(async () => {
    grok.shell.closeAll();
    grok.shell.closeTable(smiles);
  });
});
