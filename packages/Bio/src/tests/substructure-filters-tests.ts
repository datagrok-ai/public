import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, before, category, test, expect, delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {awaitGrid, readDataframe} from './utils';
import {BioSubstructureFilter, HelmFilter, SeparatorFilter} from '../widgets/bio-substructure-filter';

import {_package} from '../package-test';


category('substructureFilters', async () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadLibraries(true); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadLibraries(true); // load user settings libraries
  });

  test('fasta', async () => {
    const fasta = await readDataframe('tests/filter_FASTA.csv');
    const filter = new BioSubstructureFilter();
    await grok.data.detectSemanticTypes(fasta);
    filter.attach(fasta);
    filter.bioFilter!.substructure = 'MD';
    await delay(100);
    expect(filter.dataFrame!.filter.trueCount, 3);
    expect(filter.dataFrame!.filter.get(0), true);
    expect(filter.dataFrame!.filter.get(3), true);
    expect(filter.dataFrame!.filter.get(8), true);
    expect(filter.dataFrame!.filter.get(1), false);
  });

  test('separator', async () => {
    const msa = await readDataframe('tests/filter_MSA.csv');
    const filter = new BioSubstructureFilter();
    await grok.data.detectSemanticTypes(msa);
    filter.attach(msa);
    filter.bioFilter!.substructure = 'meI';
    await delay(100);
    expect(filter.dataFrame!.filter.trueCount, 7);
    expect(filter.dataFrame!.filter.get(2), false);
    filter.bioFilter!.substructure = '/meI';
    await delay(100);
    expect(filter.dataFrame!.filter.trueCount, 0);
    filter.bioFilter!.substructure = 'meI-hHis';
    (filter.bioFilter! as SeparatorFilter).separatorInput.value = '-';
    await delay(100);
    expect(filter.dataFrame!.filter.trueCount, 7);
    expect(filter.dataFrame!.filter.get(2), false);
  });

  test('helm', async () => {
    const helm = await readDataframe('tests/filter_HELM.csv');
    const helmTableView = grok.shell.addTableView(helm);
    const filter = new BioSubstructureFilter();
    await grok.data.detectSemanticTypes(helm);

    _package.logger.debug('Bio/substructureFilters/helm, filter attaching.');
    filter.attach(helm);
    _package.logger.debug('Bio/substructureFilters/helm, filter attached.');

    _package.logger.debug('Bio/substructureFilters/helm, filter 1 changing.');
    (filter.bioFilter! as HelmFilter).helmSubstructure = 'PEPTIDE1{C}$$$$V2.0';

    _package.logger.debug('Bio/substructureFilters/helm, filter 1 change awaiting.');
    await testEvent(helm.onRowsFiltered, () => {},
      () => { filter.bioFilter!.onChanged.next(); }, 20000);
    _package.logger.debug('Bio/substructureFilters/helm, filter 1 changed.');
    expect(filter.dataFrame!.filter.trueCount, 2);
    expect(filter.dataFrame!.filter.get(0), true);
    expect(filter.dataFrame!.filter.get(3), true);

    _package.logger.debug('Bio/substructureFilters/helm, filter 2 changing.');
    (filter.bioFilter! as HelmFilter).helmSubstructure = 'PEPTIDE1{A.C}$$$$V2.0';
    _package.logger.debug('Bio/substructureFilters/helm, filter 2 change awaiting.');
    await testEvent(helm.onRowsFiltered, () => {},
      () => { filter.bioFilter!.onChanged.next(); }, 20000);
    await awaitGrid(helmTableView.grid);
    _package.logger.debug('Bio/substructureFilters/helm, filter 2 changed.');
    expect(filter.dataFrame!.filter.trueCount, 1);
    expect(filter.dataFrame!.filter.get(3), true);
  }, {timeout: 30000});
});
