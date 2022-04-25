import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import {assessDruglikeness, drugLikenessWidget} from '../widgets/drug-likeness';
import {getIdMap, identifiersWidget} from '../widgets/identifiers';
import {getPanelElements, molfileWidget} from '../widgets/molfile';
import {getPropertiesMap, propertiesWidget} from '../widgets/properties';
import {getStructuralAlerts, structuralAlertsWidget} from '../widgets/structural-alerts';
import {structure2dWidget} from '../widgets/structure2d';
import {structure3dWidget} from '../widgets/structure3d';
import {getRisks, toxicityWidget} from '../widgets/toxicity';
import * as utils from './utils';
import $ from 'cash-dom';

category('Chem: Widgets', async () => {
  const molStr = 'CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO';

  //Test smiles, mol2000, mol3000; Check results for testing example
  test('drug-likeness', async () => {
    const dl = assessDruglikeness(molStr);
    const expectedDescription = await utils.loadFileAsText('tests/drug-likeness.json');

    expect(dl[0], 7.210692227408717);
    expect(JSON.stringify(dl[1]), expectedDescription);

    drugLikenessWidget(molStr);
  });

  //Test smiles, mol2000, mol3000; Check results for testing example
  test('identifiers', async () => {
    //RdKit Module is not initialized
    const idMap = await getIdMap(molStr);
    const expectedIdMap = await utils.loadFileAsText('tests/identifiers.json');
    expect(JSON.stringify(idMap), expectedIdMap);

    await identifiersWidget(molStr);
  });

  // Test smiles, mol2000, mol3000; Compare with existing mol string; Test copy feature
  test('molfile', async () => {
    const panelElements = getPanelElements(molStr);
    const expectedStr = (await utils.loadFileAsText('tests/molfile.sdf')).trim();
    expect(($(panelElements[2].input).val() as string).trim(), expectedStr);

    //NotAllowedError: Document is not focused.
    $(panelElements[0]).trigger('click');
    expect((await navigator.clipboard.readText()).trim(), expectedStr);

    molfileWidget(molStr);
  });

  //Test smiles, mol2000, mol3000;Compare the calculated values
  test('properties', async () => {
    const propertiesMap = getPropertiesMap(molStr);
    const expectedPropertiesMap = await utils.loadFileAsText('tests/properties.json');
    expect(JSON.stringify(propertiesMap), expectedPropertiesMap);

    propertiesWidget(molStr);
  });

  //Test smiles, mol2000, mol3000; Compare the found substructures; Visual test required
  test('structural-alerts', async () => {
    //Bad state: Cannot use origin without a scheme: undefinedfiles/alert-collection.csv
    const structuralAlerts = await getStructuralAlerts(molStr);
    const expectedSA = [1029, 1229];
    expect(structuralAlerts, expectedSA);

    await structuralAlertsWidget(molStr);
  });

  //Test smiles, mol2000, mol3000; Check if image is returned; Visual test required
  test('structure-2d', async () => {
    await grok.functions.call('structure2d', {smiles: molStr});
  });

  //Test smiles, mol2000, mol3000; Visual test required
  test('structure-3d', async () => {
    //Errors calling structure3d: molecule: Value not defined.
    await grok.functions.call('structure3d', {smiles: molStr});
  });

  //Test smiles, mol2000, mol3000; Check results for testing example
  test('toxicity', async () => {
    const risks = getRisks(molStr);
    const expectedRisks = {
      'Mutagenicity': 'None',
      'Tumorigenicity': 'None',
      'Irritating effects': 'Low',
      'Reproductive effects': 'High',
    };
    expect(JSON.stringify(risks), JSON.stringify(expectedRisks));

    toxicityWidget(molStr);
  });

  test('substructure-filter-manual', async () => {
    const df = grok.data.demo.molecules(1000);
    await grok.data.detectSemanticTypes(df);
    // previously: let filter = await grok.functions.call("Chem:substructureFilter");
    //@ts-ignore
    const filter = chem.substructureFilter();
    filter.attach(df);
    grok.shell.addTableView(df);
    const colChoice = ui.columnInput('Column', filter.dataFrame, filter.column, (col: DG.Column) => {
      filter.column = col;
      filter.dataFrame.filter.setAll(true, false);
      filter.dataFrame.rows.requestFilter();
    });
    ui.dialog({title: 'Chem Filter'})
      .add(colChoice)
      .add(filter.root)
      .show();
  });

  test('substructure-filter-panel', async () => {
    const df = DG.DataFrame.fromColumns([grok.data.demo.molecules(1000).columns[0]]);
    const view = grok.shell.addTableView(df);
    view.filters();
    await delay(100);
    (document.getElementsByClassName(
      'panel-titlebar disable-selection panel-titlebar-tabhost')[0].childNodes[0] as HTMLElement).click();
    await delay(100);
    const menuItem = document.getElementsByClassName(
      'd4-menu-item-container d4-vert-menu d4-menu-popup')[0].childNodes[3];
    menuItem.dispatchEvent(new MouseEvent('mouseenter')); await delay(100);
    menuItem.dispatchEvent(new MouseEvent('mousemove')); await delay(100);
    const submenuItem = menuItem.childNodes[0].childNodes[0].childNodes[1];
    submenuItem.dispatchEvent(new MouseEvent('mouseenter')); await delay(100);
    submenuItem.dispatchEvent(new MouseEvent('mousedown')); await delay(100);
    // https://developer.chrome.com/docs/devtools/console/utilities/#getEventListeners-function
    const dialogContents = document.getElementsByClassName('d4-dialog-contents')[0];
    const allButton = dialogContents.childNodes[2].childNodes[1].childNodes[0];
    (allButton as HTMLElement).click();
    const dialogFooter = document.getElementsByClassName('d4-dialog-footer')[0];
    const okButton = dialogFooter.childNodes[0].childNodes[0];
    (okButton as HTMLElement).click(); await delay(100);
    const smilesInput = document.getElementsByClassName('grok-sketcher-input ui-div')[0].childNodes[0];
    (smilesInput as HTMLInputElement).value = 'c1ccccc1'; await delay(100);
    smilesInput.dispatchEvent(new Event('focus')); await delay(100);
    // Only this combination of parameters worked:
    // https://tutorial.eyehunts.com/js/how-to-press-enter-key-programmatically-in-javascript-example-code/
    smilesInput.dispatchEvent(new KeyboardEvent('keydown', {
      altKey: false, bubbles: true, cancelable: true,
      charCode: 0, code: 'Enter', composed: true, ctrlKey: false,
      detail: 0, isComposing: false, key: 'Enter', keyCode: 13, location: 0,
      metaKey: false, repeat: false, shiftKey: false}));
    await delay(100);
    expect(df.filter.trueCount, 700);
  });
});
