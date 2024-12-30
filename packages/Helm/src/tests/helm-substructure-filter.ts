/* eslint-disable max-lines */
/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, delay, testEvent} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {awaitGrid, readDataframe} from './utils';

import {BioFilterProps, IBioFilter}
  from '@datagrok-libraries/bio/src/substructure-filter/bio-substructure-filter-types';

import {_package} from '../package-test';

category('bio-substructure-filters', async () => {
  let seqHelper: ISeqHelper;
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });


  test('helm-dialog', async () => {
    const logPrefix = 'Bio tests: substructureFilters/helm-dialog';
    const df = await readDataframe('tests/filter_HELM.csv');
    const view = grok.shell.addTableView(df);
    await grok.data.detectSemanticTypes(df);
    await df.meta.detectSemanticTypes();

    _package.logger.debug(`${logPrefix}, filter attaching.`);
    const filter = await grok.functions.call('Bio:bioSubstructureFilterTest');
    filter.attach(df);
    const dlg = ui.dialog('Test filters').add(filter.root).show(); // to waitForElementInDom
    await filter.awaitRendered();
    try {
      const bf = filter.bioFilter as IBioFilter;
      expect(filter.bioFilter !== null, true, 'bioFilter is not created');

      // filter 1
      _package.logger.debug(`${logPrefix}, filter 1 change awaiting...`);
      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf.props = new BioFilterProps('PEPTIDE1{A.C}$$$$V2.0', undefined, _package.logger);
      }, 20000);
      _package.logger.debug(`${logPrefix}, filter 1 changed.`);
      expect(filter.dataFrame!.filter.trueCount, 1);
      expect(filter.dataFrame!.filter.toBinaryString(), '0001');

      // filter 2
      _package.logger.debug(`${logPrefix}, filter 2 change awaiting...`);
      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf.props = new BioFilterProps('PEPTIDE1{C}$$$$V2.0', undefined, _package.logger);
      }, 20000);
      setTimeout(() => view.grid.invalidate(), 500);
      await awaitGrid(view.grid);
      await delay(1000);
      _package.logger.debug(`${logPrefix}, filter 2 changed.`);
      expect(filter.dataFrame!.filter.trueCount, 2);
      expect(filter.dataFrame!.filter.toBinaryString(), '1001');
    } finally {
      dlg.close();
    }
    await filter.awaitRendered();
    await delay(3000); //TODO: await for grid.onLookChanged
  }, {});


  // Generates unhandled exception accessing isFiltering before bioFilter created
  test('helm-view', async () => {
    const logPrefix = 'Bio tests: substructureFilters/helm-view';
    const df = await readDataframe('tests/filter_HELM.csv');
    const col = df.getCol('HELM string');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);

    const fg = view.getFiltersGroup();

    // await awaitCheck(() => fg.filters.length == 1, 'await filters.length == 1', 1000);
    // const filter = fg.filters.filter((f) => f.columnName == col.name)[0] as BioSubstructureFilter;
    await awaitGrid(view.grid);
  });


  test('sync-helm', async () => {
    const df = await _package.files.readCsv('tests/filter_HELM.csv');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);

    const fSubStr: string = 'PEPTIDE1{A.C}$$$$V2.0';
    const fTrueCount: number = 1;

    const f1 = await createFilter('HELM string', df);
    const f2 = await createFilter('HELM string', df);
    const dlg = ui.dialog('Test filters').add(f1.root).add(f2.root).show(); // to waitForElementInDom
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
    try {
      expect(!!f1.bioFilter, true);
      expect(!!f2.bioFilter, true);
      expect(f1.bioFilter!.type, 'HelmBioFilter');
      expect(f2.bioFilter!.type, 'HelmBioFilter');
      const bf1 = f1.bioFilter as IBioFilter;
      const bf2 = f2.bioFilter as IBioFilter;

      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf1.props = new BioFilterProps(fSubStr, undefined, _package.logger);
      }, 60000, 'await onRowsFiltered'); // wait to load monomers
      await awaitGrid(view.grid);
      //debugger;

      _package.logger.debug('Bio tests: substructureFilters/sync-helm, before changed event');
      await delay(f1.debounceTime * 2);
      _package.logger.debug('Bio tests: substructureFilters/sync-helm, after changed event');
      expect(df.filter.trueCount, fTrueCount);

      await f1.awaitRendered();
      expect((bf2.props as BioFilterProps).substructure, fSubStr);
    } finally {
      f1.detach();
      f2.detach();
      dlg.close();
    }
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
    await awaitGrid(view.grid);
    await delay(3000); //TODO: await for grid.onLookChanged
  });


  async function createFilter(colName: string, df: DG.DataFrame): Promise<any> {
    if (!df.columns.names().includes(colName)) {
      throw new Error(`The column '${colName}' not found. ` +
        `Available in data frame are ${JSON.stringify(df.columns.names())}`);
    }

    const filter = await grok.functions.call('Bio:bioSubstructureFilterTest');
    filter.attach(df);
    filter.applyState({columnName: colName});
    filter.column = df.col(colName);
    filter.columnName = colName;
    //filter.tableName = df.name;
    return filter;
  };
} );

