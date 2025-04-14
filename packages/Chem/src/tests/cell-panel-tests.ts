import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, test, expect, expectFloat, before, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/utils/src/test-container-utils';
import {assessDruglikeness, drugLikenessWidget} from '../widgets/drug-likeness';
import {getIdentifiersSingle} from '../widgets/identifiers';
// import {getPanelElements, molfileWidget} from '../widgets/molfile';
import {getPropertiesMap, propertiesWidget} from '../widgets/properties';
import {getStructuralAlerts, structuralAlertsWidget} from '../widgets/structural-alerts';
import {getRisks, toxicityWidget} from '../widgets/toxicity';
// import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import * as utils from './utils';
// import $ from 'cash-dom';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import * as CONST from './const';
import {structure2dWidget} from '../widgets/structure2d';
import {structure3dWidget} from '../widgets/structure3d';
import {molV2000, molV3000} from './utils';
import {EMPTY_MOLECULE_MESSAGE} from '../constants';
import {checkPackage} from '../utils/elemental-analysis-utils';
import { identifiers } from '../package';
import { getDescriptorsPy } from '../scripts-api';

const identifiersVals: {[key: string]: string} = {
  'Smiles': 'c1ccc2c(c1)CCNC2',
  'Inchi': 'InChI=1S/C9H11N/c1-2-4-9-7-10-6-5-8(9)3-1/h1-4,10H,5-7H2',
  'Inchi key': 'UWYZHKAOTLEWKK-UHFFFAOYSA-N'};

const additionalIdentifiers: {[key: string]: string} = {
  'Chembl': 'CHEMBL14346', 'PubChem': '7046', 'bindingdb': '13016',
};

category('cell panel', async () => {
  const molStr = 'CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO';
  const molFormats = [CONST.SMILES, CONST.MOL2000, CONST.MOL3000, CONST.SMARTS, CONST.EMPTY];
  const molFormats_: {[key: string]: string} = {smiles: CONST.SMILES, molV2000: CONST.MOL2000,
    molV3000: CONST.MOL3000, smarts: CONST.SMARTS, empty: CONST.EMPTY};
  const molsExt: {[key: string]: string} = {smiles: CONST.SMILES_EXTENDED, molV2000: CONST.MOL2000_EXTENDED,
      molV3000: CONST.MOL3000_EXTENDED};

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('drug-likeness', async () => {
    
    const expectedDescription = await utils.loadFileAsText('tests/drug-likeness.json');

    for (const format of Object.keys(molFormats_)) {
      if (format !== 'smarts' && format !== 'empty') {
        const dl = assessDruglikeness(molFormats_[format]);
        expectFloat(dl[0], -0.12433218410831226, undefined, `Error for ${format} format`);
        expect(JSON.stringify(dl[1]), expectedDescription, `Error for ${format} format`);
      }        
    }

    //check that widget doesn't throw errors
    for (const mol of molFormats)
      drugLikenessWidget(DG.SemanticValue.fromValueType(mol, DG.SEMTYPE.MOLECULE));
  });

  test('identifiers', async () => {
    for (const format of Object.keys(molFormats_)) {
      if (format !== 'smarts' && format !== 'empty') {
        const res: any = await getIdentifiersSingle(molFormats_[format]);
        for (const key of Object.keys(identifiersVals))
          expect(res[key], identifiersVals[key]);
        if (checkPackage('ChemblApi', 'getCompoundsIds')) {
          Object.keys(additionalIdentifiers).forEach((it) => {
            expect(res[it].innerText, additionalIdentifiers[it]);
          });
        }
        if (checkPackage('PubchemApi', 'GetIupacName'))
          expect(Object.keys(res).includes('Name'), true);
      }
    }
    //check that widget doesn't throw errors
    for (const mol of molFormats)
      identifiers(mol);
  }, { timeout: 90000 });

  test('properties', async () => {
    const expectedRes: {[key: string]: any} = {
      'MW': 133.089149,
      'HBA': 1,
      'HBD': 1,
      'LogP': 1.1253999918699265,
      'LogS': -1.8409999758005142,
      'PSA': 12.029999732971191,
      'Rotatable bonds': 0,
      'Stereo centers': 0,
      //'Molecule charge': 0,
    };
    for (const format of Object.keys(molFormats_)) {
      if (format !== 'smarts' && format !== 'empty') {
        const propertiesMap = getPropertiesMap(DG.SemanticValue.fromValueType(molFormats_[format], DG.SEMTYPE.MOLECULE));
        Object.keys(expectedRes).forEach((prop) => expect(propertiesMap[prop].value, expectedRes[prop], `Error for ${format}, ${prop} property,
          expected ${expectedRes[prop]}, got ${propertiesMap[prop].value}`));
      }
    }

    //check that widget doesn't throw errors
    for (const mol of molFormats)
      propertiesWidget(DG.SemanticValue.fromValueType(mol, DG.SEMTYPE.MOLECULE));
  });

  test('structural-alerts', async () => {
    const expectedStructuralAlerts = [1029, 1229];
    for (const format of Object.keys(molsExt)) {
      if (format !== 'smarts' && format !== 'empty') {
        const structuralAlerts = await getStructuralAlerts(molsExt[format]);
        expect(structuralAlerts.length, expectedStructuralAlerts.length);
        for (const expectedSA of expectedStructuralAlerts)
          expect(structuralAlerts.includes(expectedSA), true);
      }
    }
    //check that widget doesn't throw errors
    for (const mol of molFormats)
      await structuralAlertsWidget(mol);
  }, { timeout: 60000 });

  for (const k of Object.keys(molFormats_)) {
    test('structure2d-widget.' + k, async () => {
      structure2dWidget(molFormats_[k]);
    });
  }

  //TODO: Check if image is returned; Visual test required
  test('structure3d-widget', async () => {
    for (const mol of molFormats)
      await structure3dWidget(mol);
  });


  test('toxicity', async () => {
    const expectedRisks = {
      'Mutagenicity': 'None',
      'Tumorigenicity': 'None',
      'Irritating effects': 'Low',
      'Reproductive effects': 'High',
    };
    for (const format of Object.keys(molsExt)) {
      if (format !== 'smarts' && format !== 'empty') {
        const risks = getRisks(molsExt[format]);
        expect(JSON.stringify(risks), JSON.stringify(expectedRisks));
      }        
    }

    //check that widget doesn't throw errors
    for (const mol of molFormats)
      toxicityWidget(DG.SemanticValue.fromValueType(mol, DG.SEMTYPE.MOLECULE));
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

  test('gasteiger-partial-charges.smiles', async () => {
    const parameters = {mol: molStr, contours: 10};
    await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', parameters);
  }, {stressTest: true});

  test('gasteiger-partial-charges.molV2000', async () => {
    const parameters = {mol: molV2000, contours: 10};
    await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', parameters);
  }, {stressTest: true});

  test('gasteiger-partial-charges.molV3000', async () => {
    const parameters = {mol: molV3000, contours: 10};
    await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', parameters);
  }, {stressTest: true});

  //TODO: Compare the calculated values
  test('chem-descriptors', async () => {
    await ensureContainerRunning('name = "chem-chem"', utils.CONTAINER_TIMEOUT);
    const selesctedDesc = ["FractionCSP3", "HeavyAtomCount", "NHOHCount"];
    let jupyterRunning = false;
    grok.functions.call('Chem:TestPythonRunning', {x: 1, y: 2}).then(() => {
      console.log('*********** test python script completed');
      jupyterRunning = true; 
    });
    //check that JKG is running
    await awaitCheck(() => jupyterRunning === true, `JKG env has not been created in 2 minutes`, 120000);
    console.log('*********** chem-descriptors: started chem descriptors python script');
    const res: DG.DataFrame = await getDescriptorsPy(
      'smiles', DG.DataFrame.fromCsv(`smiles\n${molStr}`), 'selected',
      DG.DataFrame.fromColumns([DG.Column.fromList('string', 'selected', selesctedDesc)]),
    );
    console.log('*********** chem-descriptors: finished chem descriptors python script');
    expect(res.columns.names().length, 3);
    expect(selesctedDesc.every((it => res.columns.names().includes(it))), true);
    expect((res.get('FractionCSP3', 0) as number).toFixed(4), '0.6000');
    expect(res.get('HeavyAtomCount', 0), 25);
    expect(res.get('NHOHCount', 0), 1);

    //check that widget doesn't throw errors
    for (const mol of molFormats)
      await grok.functions.call('Chem:descriptorsWidget', {smiles: mol});
  }, {timeout: 60000 + utils.CONTAINER_TIMEOUT, stressTest: true});
});
