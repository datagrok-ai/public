import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, before, after, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {PharmFamilySet, runPharmacophoreDetection} from '../panels/pharmacophore-features';
import {getPharmacophoreFeatures, pharmacophoreFeaturesWidget} from '../widgets/pharmacophore-features';
import * as CONST from './const';


category('pharmacophore features', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  // Test single molecule detection — aspirin has acceptors, donors, negative, aromatic, hydrophobic
  test('single molecule: aspirin', async () => {
    const aspirin = 'CC(=O)Oc1ccccc1C(=O)O';
    const features = await getPharmacophoreFeatures(aspirin);
    expect(features.length > 0, true, 'aspirin should have pharmacophore features');
    const families = new Set(features.map((f) => f.family));
    expect(families.has('A'), true, 'aspirin should have acceptor features');
    expect(families.has('a'), true, 'aspirin should have aromatic features');
  });

  // Test donor detection — ethanolamine has O-H and N-H donors
  test('single molecule: donors', async () => {
    const ethanolamine = 'NCCO';
    const features = await getPharmacophoreFeatures(ethanolamine);
    const donors = features.filter((f) => f.family === 'D');
    expect(donors.length > 0, true, 'ethanolamine should have donor features');
  });

  // Test acceptor detection — benzonitrile has nitrile N acceptor
  test('single molecule: nitrile acceptor', async () => {
    const benzonitrile = 'N#Cc1ccccc1';
    const features = await getPharmacophoreFeatures(benzonitrile);
    const nitriles = features.filter((f) => f.featureName === 'Nitrile nitrogen');
    expect(nitriles.length > 0, true, 'benzonitrile should have nitrile nitrogen acceptor');
  });

  // Test positive ionizable — benzylamine has a primary amine
  test('single molecule: primary amine', async () => {
    const benzylamine = 'NCc1ccccc1';
    const features = await getPharmacophoreFeatures(benzylamine);
    const primaryAmines = features.filter((f) => f.featureName === 'Primary amine');
    expect(primaryAmines.length > 0, true, 'benzylamine should have primary amine feature');
  });

  // Test imidazole detection — histamine
  test('single molecule: imidazole', async () => {
    const histamine = 'NCCc1cnc[nH]1';
    const features = await getPharmacophoreFeatures(histamine);
    const imidazoles = features.filter((f) => f.featureName === 'Imidazole');
    expect(imidazoles.length > 0, true, 'histamine should have imidazole feature');
  });

  // Test negative ionizable — benzoic acid has carboxylate
  test('single molecule: negative ionizable', async () => {
    const benzoicAcid = 'OC(=O)c1ccccc1';
    const features = await getPharmacophoreFeatures(benzoicAcid);
    const negatives = features.filter((f) => f.family === 'N');
    expect(negatives.length > 0, true, 'benzoic acid should have negative ionizable features');
  });

  // Test halogen bond — chlorobenzene has halogen bond donor
  test('single molecule: halogen bond', async () => {
    const chlorobenzene = 'Clc1ccccc1';
    const features = await getPharmacophoreFeatures(chlorobenzene);
    const halogens = features.filter((f) => f.family === 'X');
    expect(halogens.length > 0, true, 'chlorobenzene should have halogen bond feature');
  });

  // Test hydrophobic — cyclohexane
  test('single molecule: hydrophobic', async () => {
    const cyclohexane = 'C1CCCCC1';
    const features = await getPharmacophoreFeatures(cyclohexane);
    const hydrophobics = features.filter((f) => f.family === 'H');
    expect(hydrophobics.length > 0, true, 'cyclohexane should have hydrophobic features');
  });

  // Test multiple molecule formats
  test('molecule formats', async () => {
    const molFormats = [CONST.SMILES, CONST.MOL2000, CONST.MOL3000];
    for (const mol of molFormats) {
      const features = await getPharmacophoreFeatures(mol);
      expect(features.length > 0, true, 'molecule should have pharmacophore features');
    }
  });

  // Test widget doesn't throw errors
  test('widget: no errors', async () => {
    const molFormats = [CONST.SMILES, CONST.MOL2000, CONST.MOL3000];
    for (const mol of molFormats)
      await pharmacophoreFeaturesWidget(mol);
  }, {timeout: 60000});

  // Test widget with empty/malformed molecules
  test('widget: empty molecule', async () => {
    const widget = await pharmacophoreFeaturesWidget('');
    expect(widget != null, true, 'widget should handle empty molecule');
  });

  // Test batch detection
  test('batch detection', async () => {
    const rdkitService = await chemCommonRdKit.getRdKitService();
    const featuresDf = DG.DataFrame.fromCsv(await _package.files.readAsText('pharmacophore-features.csv'));
    const df = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(df);
    const smilesCol = df.getCol('smiles');

    const familySet: PharmFamilySet = {
      'Donor': true, 'Acceptor': true, 'Hydrophobic': true,
      'Aromatic': true, 'Positive': false, 'Negative': false, 'Halogen Bond': false,
    };
    const resultDf = await runPharmacophoreDetection(smilesCol, familySet, featuresDf, rdkitService);
    expect(resultDf.columns.length, 4, 'should have 4 family columns (Donor, Acceptor, Hydrophobic, Aromatic)');
    expect(resultDf.rowCount, df.rowCount, 'result should have same row count as input');
  }, {timeout: 120000});

  // Test batch detection with all families enabled
  test('batch detection: all families', async () => {
    const rdkitService = await chemCommonRdKit.getRdKitService();
    const featuresDf = DG.DataFrame.fromCsv(await _package.files.readAsText('pharmacophore-features.csv'));
    const df = grok.data.demo.molecules(20);
    await grok.data.detectSemanticTypes(df);
    const smilesCol = df.getCol('smiles');

    const familySet: PharmFamilySet = {
      'Donor': true, 'Acceptor': true, 'Hydrophobic': true,
      'Aromatic': true, 'Positive': true, 'Negative': true, 'Halogen Bond': true,
    };
    const resultDf = await runPharmacophoreDetection(smilesCol, familySet, featuresDf, rdkitService);
    expect(resultDf.columns.length, 7, 'should have 7 family columns');
  }, {timeout: 120000});

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
});
