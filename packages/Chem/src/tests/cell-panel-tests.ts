import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, test, expect, expectFloat, before, waitUntil} from '@datagrok-libraries/test/src/test';
import {ensureContainerRunning} from '@datagrok-libraries/test/src/test-container-utils';
import {QED_PROPERTIES, assessDruglikeness, cruxQed, drugLikenessWidget} from '../widgets/drug-likeness';
import {getIdentifiersSingle} from '../widgets/identifiers';
// import {getPanelElements, molfileWidget} from '../widgets/molfile';
import {CRUX_PROPERTY_METHODS, getPropertiesMap, propertiesWidget} from '../widgets/properties';
import {ChemicalEngine} from '../engine/chem-engine';
import {inEngine} from './engine-utils';
import {getStructuralAlerts, structuralAlertsWidget} from '../widgets/structural-alerts';
import {cruxToxicityAlerts, getRisks, toxicityWidget} from '../widgets/toxicity';
// import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import * as utils from './utils';
// import $ from 'cash-dom';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import * as CONST from './const';
import {structure2dWidget} from '../widgets/structure2d';
import {cruxEmbed3d, embed3dOptions, isEmbedTimeout, structure3dWidget} from '../widgets/structure3d';
import {_convertMolNotation} from '../utils/convert-notation-utils';
import {molV2000, molV3000} from './utils';
import {checkPackage} from '../utils/elemental-analysis-utils';
import {PackageFunctions} from '../package';
import {calculateDescriptors} from '../docker/api';

const identifiersVals: {[key: string]: string} = {
  'Smiles': 'c1ccc2ccccc2c1',
  'Inchi': 'InChI=1S/C10H8/c1-2-6-10-8-4-3-7-9(10)5-1/h1-8H',
  'Inchi key': 'UFWIBTONFRDIAS-UHFFFAOYSA-N'};

const additionalIdentifiers: { [key: string]: string } = {
  'Chembl': 'CHEMBL16293',
  'PubChem': '931',
  'bindingdb': '50159249',
  'pdb': 'NPY',
  'chebi': 'CHEBI:16482',
  'fdasrs': '2166IN72UN',
  'hmdb': 'HMDB0029751',
  'brenda': '1830',
  'surechembl': '8953',
};

category('cell panel', async () => {
  const molStr = 'CC(C)Cc1ccc(cc1)C(C)C(=O)N2CCCC2C(=O)OCCO';
  const molFormats = [CONST.SMILES, CONST.MOL2000, CONST.MOL3000, CONST.SMARTS, CONST.EMPTY];
  const molFormats_: {[key: string]: string} = {smiles: CONST.SMILES, molV2000: CONST.MOL2000,
    molV3000: CONST.MOL3000, smarts: CONST.SMARTS, empty: CONST.EMPTY};
  const molFormatsForIdentifiers_: {[key: string]: string} = {smiles: CONST.SMILES_2, molV2000: CONST.MOL2000_2,
    molV3000: CONST.MOL3000_2};
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

  test('drug-likeness: Crux mode, QED', async () => {
    // RDKit's QED.qed and QED.properties (the pinned RDKit's values, crux-core's tools/chem-wasm/
    // chem_qed_alerts_oracle.py), not OpenChemLib's score: the product owner, 2026-10-22
    const expected: {[key: string]: number} = {'MW': 133.194, 'ALOGP': 1.3323, 'HBA': 1, 'HBD': 1, 'PSA': 12.03,
      'ROTB': 0, 'AROM': 1, 'ALERTS': 0};
    await inEngine(ChemicalEngine.Crux, async () => {
      for (const format of ['smiles', 'molV2000', 'molV3000']) {
        const qed = cruxQed(molFormats_[format]);
        expect(qed !== null, true, `${format}: read`);
        expectFloat(qed!['QED'], 0.5627155967397938, 1e-12, `${format}: QED`);
        for (const [key, value] of Object.entries(expected))
          expect(qed![key], value, `${format}: ${key}`);
      }
      const widget = drugLikenessWidget(DG.SemanticValue.fromValueType(CONST.SMILES, DG.SEMTYPE.MOLECULE));
      expect(widget.root.textContent!.includes('QED: 0.56'), true, 'the widget shows QED');
      const labels = Array.from(widget.root.querySelectorAll('tr')).map((tr) => tr.cells[0]?.textContent);
      expect(labels.join(', '), Object.values(QED_PROPERTIES).join(', '), 'and the properties it weighs');
    });
  });

  test('identifiers', async () => {
    for (const format of Object.keys(molFormatsForIdentifiers_)) {
      if (format !== 'smarts' && format !== 'empty') {
        const res: any = await getIdentifiersSingle(molFormatsForIdentifiers_[format]);
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
      PackageFunctions.identifiers(mol);
  }, {timeout: 90000});

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
    await inEngine(ChemicalEngine.RDKit, async () => { // OpenChemLib's values (Crux mode's: the next test)
      for (const format of Object.keys(molFormats_)) {
        if (format !== 'smarts' && format !== 'empty') {
          const propertiesMap =
            getPropertiesMap(DG.SemanticValue.fromValueType(molFormats_[format], DG.SEMTYPE.MOLECULE));
          Object.keys(expectedRes).forEach((prop) => expect(propertiesMap[prop].value, expectedRes[prop],
            `Error for ${format}, ${prop} property, expected ${expectedRes[prop]}, got ${propertiesMap[prop].value}`));
        }
      }

      //check that widget doesn't throw errors
      for (const mol of molFormats)
        propertiesWidget(DG.SemanticValue.fromValueType(mol, DG.SEMTYPE.MOLECULE));
    });
  });

  test('properties: Crux mode, RDKit\'s definitions', async () => {
    // crux's chemProperties: RDKit's values (MolWt, Lipinski HBA, NumHDonors, Crippen, ESOL, TPSA with S and P, ...;
    // the expected values from the pinned RDKit, crux-core's tools/chem-wasm/chem_props_oracle.py), not OpenChemLib's
    const expected: {[key: string]: any} = {'MW': 133.194, 'HBA': 1, 'HBD': 1, 'LogP': 1.3323,
      'LogS': -1.8636003565440682, 'PSA': 12.03, 'Rotatable bonds': 0, 'Stereo centers': 0, 'Molecule charge': 0,
      'Molecular formula': 'C9H11N'};
    await inEngine(ChemicalEngine.Crux, async () => {
      for (const format of ['smiles', 'molV2000', 'molV3000']) {
        const map = getPropertiesMap(DG.SemanticValue.fromValueType(molFormats_[format], DG.SEMTYPE.MOLECULE));
        for (const [name, value] of Object.entries(expected)) {
          expect(map[name].value, value, `${format}: ${name}`);
          expect(map[name].method, CRUX_PROPERTY_METHODS[name], `${format}: ${name}'s method`);
        }
      }
      // the widget names each property's method and reference
      const widget = propertiesWidget(DG.SemanticValue.fromValueType(CONST.SMILES, DG.SEMTYPE.MOLECULE));
      const names = Array.from(widget.root.querySelectorAll('tr')).map((tr) => tr.cells[0]?.textContent);
      expect(names.join(', '), Object.keys(expected).join(', '));
      expect(CRUX_PROPERTY_METHODS['LogS'].includes('Delaney'), true);
    });
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
  }, {timeout: 60000});

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

  test('structure3d-widget: Crux mode, chemEmbed3d', async () => {
    // the panel's 2D MOL block embedded by crux's chemEmbed3d on a Crux worker (ETKDGv3 without Hs, crux-core decision
    // chem-props), no Python; the molecule as drawn: its atoms and bonds, 3D coordinates
    const rdKitModule = chemCommonRdKit.getRdKitModule();
    expect(await inEngine(ChemicalEngine.RDKit, () => cruxEmbed3d(CONST.SMILES)), null, 'RDKit mode: the script');
    await inEngine(ChemicalEngine.Crux, async () => {
      for (const format of ['smiles', 'molV2000', 'molV3000']) {
        const flat = _convertMolNotation(molFormats_[format], DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock,
          rdKitModule);
        const block = await cruxEmbed3d(flat);
        expect(block !== null, true, `${format}: embedded`);
        const lines = (b: string) => b.split('\n');
        const counts = (b: string) => lines(b)[3].substring(0, 6);
        expect(counts(block!), counts(flat), `${format}: the atom and bond counts`);
        const n = Number(counts(flat).substring(0, 3));
        const atoms = lines(block!).slice(4, 4 + n);
        expect(atoms.map((l) => l.substring(31, 34)).join(), lines(flat).slice(4, 4 + n).map((l) =>
          l.substring(31, 34)).join(), `${format}: the atoms`);
        expect(atoms.some((l) => Number(l.substring(20, 30)) !== 0), true, `${format}: 3D coordinates`);
      }
      expect(await cruxEmbed3d(CONST.EMPTY), null, 'no block for no atoms: the script\'s path');
      // the registered entry point runs in the package's own bundle, under the stand's engine (not this test's switch)
      if (await grok.functions.call('Chem:chemicalEngine') === ChemicalEngine.Crux) {
        expect(String(await grok.functions.call('Chem:embed3d', {molecule: CONST.SMILES})).includes('M  END'), true,
          'Chem:embed3d, the entry point other packages call');
      }
      await structure3dWidget(CONST.SMILES);
    });
  });

  test('structure3d-widget: Crux mode, a large macrocycle embeds within its deadline, or gives up', async () => {
    // a macrocyclic peptide the old embedding (ETKDG v1, unbounded attempts) never finished, freezing the pane
    // (the product owner, 2026-10-23); with Crux mode's options it embeds, and a deadline too short for it throws
    // and leaves the worker usable
    const peptide = 'CCN=C(NCC)NCCCC[C@H]1NC(=O)[C@H](Cc2cc(I)c(O)c(I)c2)NC(=O)[C@@H]([C@H](C)O)NC(=O)[C@H](Cc2cccnc2' +
      ')NC(=O)[C@H](Cc2ccc(Oc3ccc(O)cc3)cc2)NC(=O)c2ccccc2N(C)C(=O)C(C)=NC(=O)[C@@H](Cc2ccc(F)cc2)NC(=O)[C@H](Cc2cscn' +
      '2)NC(=O)[C@H](C(C)C)N(C)C(=O)[C@@H](Cc2ccc(OC)cc2)NC(=O)CN(C)C(=O)[C@H]2CCCN2C(=O)Cn2ccc3ccc(cc3c2=O)N2CCN(CC2' +
      ')C(=O)[C@H]2C[C@@H](O)CN2C(=O)[C@H]([C@@H](C)O)N(C)C1=O';
    await inEngine(ChemicalEngine.Crux, async () => {
      let started = performance.now();
      const block = await cruxEmbed3d(peptide);
      expect(block !== null && block.includes('M  END'), true, 'the peptide embeds');
      expect(performance.now() - started < embed3dOptions().timeoutMs, true, 'within the deadline');
      started = performance.now();
      let error: unknown = null;
      try {
        // the old path (ETKDG v1, the legacy minimization, no random coordinates), which never embedded it
        await cruxEmbed3d(peptide, {etkdg: 'v1', legacy: true, randomCoords: false, maxIterations: 0, timeoutMs: 300});
      } catch (e) {
        error = e;
      }
      expect(isEmbedTimeout(error), true, `a 300 ms deadline throws the timeout (got ${error})`);
      expect(performance.now() - started < 300 + 1500, true, 'soon after the deadline');
      expect((await cruxEmbed3d(CONST.SMILES))?.includes('M  END'), true, 'the worker embeds the next molecule');
      const widget = await structure3dWidget(CONST.SMILES);
      expect(widget.root.textContent?.includes('No 3D coordinates'), false, 'the pane draws');
    });
  }, {timeout: 60000});


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


  test('toxicity: Crux mode, Brenk and NIH alerts', async () => {
    // RDKit's FilterCatalog Brenk and NIH sets (the pinned RDKit's hits and atoms, crux-core's tools/chem-wasm/
    // chem_qed_alerts_oracle.py), not OpenChemLib's risks: the product owner, 2026-10-22
    const thioester = 'C=CC(=O)NCCSC(=O)C';
    await inEngine(ChemicalEngine.Crux, async () => {
      expect(JSON.stringify(cruxToxicityAlerts(CONST.SMILES)), '{"Brenk":[],"NIH":[]}', 'no alert');
      const alerts = cruxToxicityAlerts(thioester)!;
      expect(JSON.stringify(alerts['Brenk']), '[{"name":"Aliphatic_long_chain","atoms":[[4,5,6,7]]},' +
        '{"name":"Michael_acceptor_1","atoms":[[0,1,2,3]]},{"name":"thioester","atoms":[[7,8,9]]}]', 'Brenk');
      expect(JSON.stringify(alerts['NIH']), '[{"name":"acrylate","atoms":[[0,1,2,3]]},' +
        '{"name":"halo_acrylate","atoms":[[0,1]]},{"name":"thioester","atoms":[[7,8,9,10]]}]', 'NIH');
      // the widget: each catalog's verdict, each alert drawn
      const widget = toxicityWidget(DG.SemanticValue.fromValueType(thioester, DG.SEMTYPE.MOLECULE));
      const verdicts = Array.from(widget.root.querySelectorAll('.chem-toxicity-widget-risk-row'))
        .map((d) => d.textContent);
      expect(verdicts.join(', '), '3 alerts, 3 alerts', 'the verdicts');
      expect(widget.root.querySelectorAll('canvas').length, 6, 'an image per alert');
      const none = toxicityWidget(DG.SemanticValue.fromValueType(CONST.SMILES, DG.SEMTYPE.MOLECULE));
      expect(Array.from(none.root.querySelectorAll('.chem-toxicity-widget-risk-row')).map((d) => d.textContent)
        .join(', '), 'None, None', 'no alert');
    });
  });


  //TODO: fix, Test smiles, mol2000, mol3000;
  // test('substructure-filter-panel', async () => {
  //   const df = DG.DataFrame.fromColumns([grok.data.demo.molecules(1000).columns.byIndex(0)]);
  //   const view = grok.shell.addTableView(df);
  //   view.filters();
  //   const titleBars = () =>
  //     document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost');
  //   await waitUntil(() => titleBars().length > 0, 'the filter panel did not open');
  //   (titleBars()[0].childNodes[0] as HTMLElement).click();
  //   const menus = () => document.getElementsByClassName('d4-menu-item-container d4-vert-menu d4-menu-popup');
  //   await waitUntil(() => menus().length > 0, 'the filter panel menu did not open');
  //   const menuItem = menus()[0].childNodes[3];
  //   menuItem.dispatchEvent(new MouseEvent('mouseenter'));
  //   menuItem.dispatchEvent(new MouseEvent('mousemove'));
  //   await waitUntil(() => menuItem.childNodes[0]?.childNodes[0]?.childNodes[1] != null, 'the submenu did not open');
  //   const submenuItem = menuItem.childNodes[0].childNodes[0].childNodes[1];
  //   submenuItem.dispatchEvent(new MouseEvent('mouseenter'));
  //   submenuItem.dispatchEvent(new MouseEvent('mousedown'));
  //   await waitUntil(() => document.getElementsByClassName('d4-dialog-contents').length > 0, 'no columns dialog');
  //   // https://developer.chrome.com/docs/devtools/console/utilities/#getEventListeners-function
  //   const dialogContents = document.getElementsByClassName('d4-dialog-contents')[0];
  //   const allButton = dialogContents.childNodes[2].childNodes[1].childNodes[0];
  //   (allButton as HTMLElement).click();
  //   const dialogFooter = document.getElementsByClassName('d4-dialog-footer')[0];
  //   const okButton = dialogFooter.childNodes[0].childNodes[0];
  //   (okButton as HTMLElement).click();
  //   const inputs = () => document.getElementsByClassName('grok-sketcher-input ui-div');
  //   await waitUntil(() => inputs().length > 0, 'the substructure filter shows no sketcher input');
  //   const smilesInput = inputs()[0].childNodes[0];
  //   (smilesInput as HTMLInputElement).value = 'c1ccccc1';
  //   smilesInput.dispatchEvent(new Event('focus'));
  //   // Only this combination of parameters worked:
  //   // https://tutorial.eyehunts.com/js/how-to-press-enter-key-programmatically-in-javascript-example-code/
  //   smilesInput.dispatchEvent(new KeyboardEvent('keydown', {
  //     altKey: false, bubbles: true, cancelable: true,
  //     charCode: 0, code: 'Enter', composed: true, ctrlKey: false,
  //     detail: 0, isComposing: false, key: 'Enter', keyCode: 13, location: 0,
  //     metaKey: false, repeat: false, shiftKey: false}));
  //   await waitUntil(() => df.filter.trueCount === 700, () => `${df.filter.trueCount} rows pass, 700 expected`);
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
    const selesctedDesc = ['FractionCSP3', 'HeavyAtomCount', 'NHOHCount'];
    const cols = await calculateDescriptors(DG.Column.fromStrings('smiles', [molStr]), selesctedDesc);
    expect(cols.length, 3);
    const colNames = cols.map((it) => it.name);
    expect(selesctedDesc.every(((it) => colNames.includes(it))), true);
    expect((cols.filter((it) => it.name === 'FractionCSP3')[0].get(0) as number).toFixed(4), '0.6000');
    expect(cols.filter((it) => it.name === 'HeavyAtomCount')[0].get(0), 25);
    expect(cols.filter((it) => it.name === 'NHOHCount')[0]!.get(0), 1);

    //check that widget doesn't throw errors
    for (const mol of molFormats)
      await grok.functions.call('Chem:descriptorsWidget', {smiles: mol});
  }, {timeout: 60000 + utils.CONTAINER_TIMEOUT, stressTest: true});
});
