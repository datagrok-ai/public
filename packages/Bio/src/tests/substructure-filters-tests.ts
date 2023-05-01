import {after, before, category, test, expect, delay} from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {readDataframe} from './utils';
import {BioSubstructureFilter, HelmFilter, SeparatorFilter} from '../widgets/bio-substructure-filter';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {LIB_DEFAULT, LIB_STORAGE_NAME} from '../utils/monomer-lib';


category('substructureFilters', async () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibrariesSettings: {};

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibrariesSettings = await grok.dapi.userDataStorage.get(LIB_STORAGE_NAME, true);

    // Test 'helm' requires default monomer library loaded
    await grok.dapi.userDataStorage.post(LIB_STORAGE_NAME, LIB_DEFAULT, true);
    await monomerLibHelper.loadLibraries(true); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await grok.dapi.userDataStorage.put(LIB_STORAGE_NAME, userLibrariesSettings, true);
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
    filter.attach(helm);

    const helmFilterChanged = new Promise((resolve, reject) => {
      helm.onFilterChanged.subscribe(async (_: any) => {
        try {
          resolve(true);
        } catch (error) {
          reject(error);
        }
      });
    });
    (filter.bioFilter! as HelmFilter).helmSubstructure = 'PEPTIDE1{C}$$$$V2.0';
    filter.bioFilter!.onChanged.next();
    await helmFilterChanged;

    //await delay(3000);
    expect(filter.dataFrame!.filter.trueCount, 2);
    expect(filter.dataFrame!.filter.get(0), true);
    expect(filter.dataFrame!.filter.get(3), true);
    (filter.bioFilter! as HelmFilter).helmSubstructure = 'PEPTIDE1{A.C}$$$$V2.0';
    filter.bioFilter!.onChanged.next();
    await delay(100);
    expect(filter.dataFrame!.filter.trueCount, 1);
    expect(filter.dataFrame!.filter.get(3), true);
    helmTableView.close();
  }, {skipReason: 'GROK-12779'});
});
