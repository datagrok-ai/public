import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, test, expect, delay, expectFloat, before} from '@datagrok-libraries/utils/src/test';
import {assessDruglikeness, drugLikenessWidget} from '../widgets/drug-likeness';
import {getIdMap, identifiersWidget} from '../widgets/identifiers';
import {getPanelElements, molfileWidget} from '../widgets/molfile';
import {getPropertiesMap, propertiesWidget} from '../widgets/properties';
import {getStructuralAlerts, structuralAlertsWidget} from '../widgets/structural-alerts';
import {getRisks, toxicityWidget} from '../widgets/toxicity';
import * as utils from './utils';
import $ from 'cash-dom';
import { _package } from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import { getDescriptorsSingle } from '../descriptors/descriptors-calculation';
import {runStructuralAlertsDetection} from '../panels/structural-alerts';
import {RDMol} from '../rdkit-api';


category('widgets', async () => {
  const molStr = 'CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO';

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  //TODO: Test mol2000, mol3000
  test('drug-likeness', async () => {
    const dl = assessDruglikeness(molStr);
    const expectedDescription = await utils.loadFileAsText('tests/drug-likeness.json');

    expectFloat(dl[0], 7.210692227408717);
    expect(JSON.stringify(dl[1]), expectedDescription);

    drugLikenessWidget(molStr);
  });

  //TODO: Test mol2000, mol3000
  test('identifiers', async () => {
    const idMap = await getIdMap(molStr);
    const expectedIdMap = await utils.loadFileAsText('tests/identifiers.json');
    expect(JSON.stringify(idMap), expectedIdMap);

    await identifiersWidget(molStr);
  });

  //TODO: Test mol2000, mol3000
  test('molfile', async () => {
    const expectedStr = (await utils.loadFileAsText('tests/molfile.sdf')).replaceAll('\r', '').trim();
    const panelElements = getPanelElements(molStr);
    const panelStr = ($(panelElements[2].input).val() as string).replaceAll('\r', '').trim();
    expect(panelStr, expectedStr);

    $(panelElements[0]).trigger('click');
    const clipboardStr = (await navigator.clipboard.readText()).replaceAll('\r', '').trim();
    expect(clipboardStr, expectedStr);

    molfileWidget(molStr);
  });

  //TODO: Test mol2000, mol3000;
  test('properties', async () => {
    const propertiesMap = getPropertiesMap(molStr);
    const expectedPropertiesMap = await utils.loadFileAsText('tests/properties.json');
    expect(JSON.stringify(propertiesMap), expectedPropertiesMap);

    propertiesWidget(molStr);
  });

  //TODO: Test mol2000, mol3000; Visual test required; Unstable
  test('structural-alerts', async () => {
    const structuralAlerts = await getStructuralAlerts(molStr);
    const expectedStructuralAlerts = [1029, 1229];
    expect(structuralAlerts.length, expectedStructuralAlerts.length);
    for (const expectedSA of expectedStructuralAlerts)
      expect(structuralAlerts.includes(expectedSA), true);

    await structuralAlertsWidget(molStr);
  });

  //TODO: Test mol2000, mol3000; Check if image is returned; Visual test required
  test('structure-2d', async () => {
    await grok.functions.call('structure2d', {smiles: molStr});
  });

  //TODO: Test mol2000, mol3000; Visual test required
  test('structure-3d', async () => {
    await grok.functions.call('structure3d', {smiles: molStr});
  });

  //TODO: Test mol2000, mol3000
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

  //TODO: Test smiles, mol2000, mol3000;
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

  //TODO: fix, Test smiles, mol2000, mol3000;
  // test('substructure-filter-panel', async () => {
  //   const df = DG.DataFrame.fromColumns([grok.data.demo.molecules(1000).columns.byIndex(0)]);
  //   const view = grok.shell.addTableView(df);
  //   view.filters();
  //   await delay(100);
  //   (document.getElementsByClassName(
  //     'panel-titlebar disable-selection panel-titlebar-tabhost')[0].childNodes[0] as HTMLElement).click();
  //   await delay(100);
  //   const menuItem = document.getElementsByClassName(
  //     'd4-menu-item-container d4-vert-menu d4-menu-popup')[0].childNodes[3];
  //   menuItem.dispatchEvent(new MouseEvent('mouseenter')); await delay(100);
  //   menuItem.dispatchEvent(new MouseEvent('mousemove')); await delay(100);
  //   const submenuItem = menuItem.childNodes[0].childNodes[0].childNodes[1];
  //   submenuItem.dispatchEvent(new MouseEvent('mouseenter')); await delay(100);
  //   submenuItem.dispatchEvent(new MouseEvent('mousedown')); await delay(100);
  //   // https://developer.chrome.com/docs/devtools/console/utilities/#getEventListeners-function
  //   const dialogContents = document.getElementsByClassName('d4-dialog-contents')[0];
  //   const allButton = dialogContents.childNodes[2].childNodes[1].childNodes[0];
  //   (allButton as HTMLElement).click();
  //   const dialogFooter = document.getElementsByClassName('d4-dialog-footer')[0];
  //   const okButton = dialogFooter.childNodes[0].childNodes[0];
  //   (okButton as HTMLElement).click(); await delay(100);
  //   const smilesInput = document.getElementsByClassName('grok-sketcher-input ui-div')[0].childNodes[0];
  //   (smilesInput as HTMLInputElement).value = 'c1ccccc1'; await delay(100);
  //   smilesInput.dispatchEvent(new Event('focus')); await delay(100);
  //   // Only this combination of parameters worked:
  //   // https://tutorial.eyehunts.com/js/how-to-press-enter-key-programmatically-in-javascript-example-code/
  //   smilesInput.dispatchEvent(new KeyboardEvent('keydown', {
  //     altKey: false, bubbles: true, cancelable: true,
  //     charCode: 0, code: 'Enter', composed: true, ctrlKey: false,
  //     detail: 0, isComposing: false, key: 'Enter', keyCode: 13, location: 0,
  //     metaKey: false, repeat: false, shiftKey: false}));
  //   await delay(100);
  //   expect(df.filter.trueCount, 700);
  // });

  //TODO: Test mol2000, mol3000;
  test('gasteiger-partion-charges', async () => {
    const parameters = {mol: molStr, contours: 10};
    await grok.functions.call('Chem:GasteigerPartialCharges', parameters);
  });

  //TODO: Test mol2000, mol3000; Compare the calculated values
  test('chem-descriptors', async () => {
    getDescriptorsSingle(molStr);
  });
});

category('widget benchmarks', () => {
  before(async () => {
    chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
    chemCommonRdKit.initRdKitModuleLocal();
  });

  test('Structural Alerts Benchmark', async () => {
    const alertsDf = DG.DataFrame.fromCsv(await _package.files.readAsText('alert-collection.csv'));
    const ruleSetCol = alertsDf.getCol('rule_set_name');
    const smartsCol = alertsDf.getCol('smarts');
    const ruleIdCol = alertsDf.getCol('rule_id');
    const rdkitModule = chemCommonRdKit.getRdKitModule();
  
    const smartsMap = new Map<string, RDMol>();
    for (let i = 0; i < alertsDf.rowCount; i++)
      smartsMap.set(ruleIdCol.get(i), rdkitModule.get_qmol(smartsCol.get(i)));

    const sarSmall = DG.DataFrame.fromCsv(await _package.files.readAsText('sar-small.csv'));
    const smilesCol = sarSmall.getCol('smiles');
    const ruleSetList = ['BMS', 'Dundee', 'Glaxo', 'Inpharmatica', 'LINT', 'MLSMR', 'PAINS', 'SureChEMBL'];

    DG.time('Structural Alerts', () => {
      runStructuralAlertsDetection(sarSmall, ruleSetList, smilesCol, ruleSetCol, ruleIdCol, smartsMap, rdkitModule);
    });
  });
});