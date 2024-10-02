import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';

import {after, before, category, test, expect, delay, testEvent, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {awaitGrid, readDataframe} from './utils';

import {
  BioSubstructureFilter, FastaBioFilter, SeparatorBioFilter, SeparatorFilterProps
} from '../widgets/bio-substructure-filter';
import {BioFilterProps} from '../widgets/bio-substructure-filter-types';
import {HelmBioFilter} from '../widgets/bio-substructure-filter-helm';

import {_package} from '../package-test';


category('bio-substructure-filters', async () => {
  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    helmHelper = await getHelmHelper(); // init Helm package
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

  test('fasta', async () => {
    const df = await readDataframe('tests/filter_FASTA.csv');
    await grok.data.detectSemanticTypes(df);

    const fSubStr: string = 'MD';
    const fTrueCount: number = 3;

    const filter = new BioSubstructureFilter();
    filter.attach(df);
    await filter.awaitRendered();
    try {
      expect(!!filter.bioFilter, true);
      expect(filter.bioFilter!.type, 'FastaBioFilter');
      const bf = filter.bioFilter as FastaBioFilter;

      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf.props = new BioFilterProps(fSubStr);
      }, 5000, 'testEvent onRowsFiltered');
      expect(filter.dataFrame!.filter.trueCount, fTrueCount);
      expect(filter.dataFrame?.filter.toBinaryString(), '10010000100000');
    } finally {
      filter.detach();
    }
    await filter.awaitRendered();
  });

  test('separator', async () => {
    const msa = await readDataframe('tests/filter_MSA.csv');
    const filter = new BioSubstructureFilter();
    await grok.data.detectSemanticTypes(msa);
    filter.attach(msa);
    await filter.awaitRendered();
    try {
      expect(!!filter.bioFilter, true);
      expect(filter.bioFilter!.type, 'SeparatorBioFilter');
      const bf = filter.bioFilter as SeparatorBioFilter;

      await testEvent(msa.onRowsFiltered, () => {}, () => {
        bf.props = new SeparatorFilterProps('meI');
      }, 5000, 'testEvent onRowsFiltered');
      expect(filter.dataFrame!.filter.trueCount, 7);
      expect(filter.dataFrame!.filter.get(2), false);

      await testEvent(msa.onRowsFiltered, () => {}, () => {
        bf.props = new SeparatorFilterProps('/meI');
      }, 5000, 'testEvent onRowsFiltered');
      expect(filter.dataFrame!.filter.trueCount, 0);

      await testEvent(msa.onRowsFiltered, () => {}, () => {
        bf.props = new SeparatorFilterProps('meI-hHis', '-');
      }, 5000, 'testEvent onRowsFiltered');
      expect(filter.dataFrame!.filter.trueCount, 7);
      expect(filter.dataFrame!.filter.get(2), false);
    } finally {
      filter.detach();
    }
    await filter.awaitRendered();
  });

  // test('helm', async () => {
  //   const df = await readDataframe('tests/filter_HELM.csv');
  //   await grok.data.detectSemanticTypes(df);
  //   const view = grok.shell.addTableView(df);
  //
  //   // // Helm filter calls waitForElementInDom
  //   // const fg = view.getFiltersGroup({createDefaultFilters: false});
  //   // _package.logger.debug('Bio tests: substructureFilters/helm, filter attaching.');
  //   // const filter = new BioSubstructureFilter();
  //   // filter.attach(df);
  //   // _package.logger.debug('Bio tests: substructureFilters/helm, filter attached.');
  //   // const fg = await df.plot.fromType(DG.VIEWER.FILTERS, {
  //   //   filters: [
  //   //     {type: 'Bio:bioSubstructureFilter', column: 'HELM string'}]
  //   // }) as DG.FilterGroup;
  //   // const fg2 = view.getFiltersGroup({createDefaultFilters: false});
  //   const fg = view.getFiltersGroup();
  //   await awaitCheck(() => fg.filters.length == 1, 'Filter panel has not been created', 3000);
  //   await awaitGrid(view.grid);
  //   const filter = fg.filters[0] as BioSubstructureFilter;
  //   // TODO: Check filter type
  //
  //   _package.logger.debug('Bio tests: substructureFilters/helm, filter 1 changing.');
  //   const bf = filter.bioFilter as HelmBioFilter;
  //   await delay(1000); // to draw grid and filters
  //   await bf.awaitRendered();
  //
  //   _package.logger.debug('Bio tests: substructureFilters/helm, filter 1 change awaiting.');
  //   await testEvent(df.onRowsFiltered, () => {}, () => {
  //     bf.props = new BioFilterProps('PEPTIDE1{C}$$$$V2.0');
  //   }, 20000);
  //   _package.logger.debug('Bio tests: substructureFilters/helm, filter 1 changed.');
  //   expect(filter.dataFrame!.filter.trueCount, 2);
  //   expect(filter.dataFrame!.filter.get(0), true);
  //   expect(filter.dataFrame!.filter.get(3), true);
  //
  //   _package.logger.debug('Bio tests: substructureFilters/helm, filter 2 changing.');
  //   await testEvent(df.onRowsFiltered, () => {}, () => {
  //     bf.props = new BioFilterProps('PEPTIDE1{A.C}$$$$V2.0');
  //   }, 20000);
  //   await awaitGrid(view.grid);
  //   _package.logger.debug('Bio tests: substructureFilters/helm, filter 2 changed.');
  //   expect(filter.dataFrame!.filter.trueCount, 1);
  //   expect(filter.dataFrame!.filter.get(3), true);
  // }, {timeout: 30000});

  test('helm-dialog', async () => {
    const logPrefix = 'Bio tests: substructureFilters/helm-dialog';
    const df = await readDataframe('tests/filter_HELM.csv');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);

    _package.logger.debug(`${logPrefix}, filter attaching.`);
    const filter = new BioSubstructureFilter();
    filter.attach(df);
    const dlg = ui.dialog('Test filters').add(filter.root).show(); // to waitForElementInDom
    await filter.awaitRendered();
    try {
      const bf = filter.bioFilter as HelmBioFilter;
      expect(filter.bioFilter !== null, true, 'bioFilter is not created');

      // filter 1
      _package.logger.debug(`${logPrefix}, filter 1 change awaiting...`);
      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf.props = new BioFilterProps('PEPTIDE1{A.C}$$$$V2.0');
      }, 20000);
      _package.logger.debug(`${logPrefix}, filter 1 changed.`);
      expect(filter.dataFrame!.filter.trueCount, 1);
      expect(filter.dataFrame!.filter.toBinaryString(), '0001');

      // filter 2
      _package.logger.debug(`${logPrefix}, filter 2 change awaiting...`);
      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf.props = new BioFilterProps('PEPTIDE1{C}$$$$V2.0');
      }, 20000);
      await awaitGrid(view.grid);
      _package.logger.debug(`${logPrefix}, filter 2 changed.`);
      expect(filter.dataFrame!.filter.trueCount, 2);
      expect(filter.dataFrame!.filter.toBinaryString(), '1001');
    } finally {
      dlg.close();
    }
    await filter.awaitRendered();
    await delay(3000); //TODO: await for grid.onLookChanged
  });

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

  test('sync-fasta', async () => {
    const df = await _package.files.readCsv('tests/filter_FASTA.csv');
    await grok.data.detectSemanticTypes(df);

    const fSubStr: string = 'MD';
    const fTrueCount: number = 3;

    const f1 = await createFilter('fasta', df);
    const f2 = await createFilter('fasta', df);
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
    try {
      expect(!!f1.bioFilter, true);
      expect(!!f2.bioFilter, true);
      expect(f1.bioFilter!.type, 'FastaBioFilter');
      expect(f2.bioFilter!.type, 'FastaBioFilter');
      const bf1 = f1.bioFilter as FastaBioFilter;
      const bf2 = f2.bioFilter as FastaBioFilter;

      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf1.props = new BioFilterProps(fSubStr);
      }, 10000, 'await onRowsFiltered');
      expect(df.filter.trueCount, fTrueCount);

      await f1.awaitRendered();
      expect(bf2.props.substructure, fSubStr);
    } finally {
      f1.detach();
      f2.detach();
    }
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
  });

  // MSA filter has the second input field for separator
  test('sync-msa', async () => {
    const df = await _package.files.readCsv('tests/filter_MSA.csv');
    await grok.data.detectSemanticTypes(df);

    const fSubStr: string = 'hHis-Aca';
    const fSepStr: string = '-';
    const fTrueCount: number = 8;

    const f1 = await createFilter('MSA', df);
    const f2 = await createFilter('MSA', df);
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
    try {
      expect(!!f1.bioFilter, true);
      expect(!!f2.bioFilter, true);
      expect(f1.bioFilter!.type, 'SeparatorBioFilter');
      expect(f2.bioFilter!.type, 'SeparatorBioFilter');
      const bf1 = f1.bioFilter as SeparatorBioFilter;
      const bf2 = f2.bioFilter as SeparatorBioFilter;

      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf1.props = new SeparatorFilterProps(fSubStr, fSepStr);
      }, 10000, 'await onRowsFiltered');
      expect(df.filter.trueCount, fTrueCount);

      expect(bf2.props.substructure, fSubStr);
      expect(bf2.props.separator, fSepStr);
    } finally {
      f1.detach();
      f2.detach();
    }
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
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
      const bf1 = f1.bioFilter as HelmBioFilter;
      const bf2 = f2.bioFilter as HelmBioFilter;

      await testEvent(df.onRowsFiltered, () => {}, () => {
        bf1.props = new BioFilterProps(fSubStr);
      }, 60000, 'await onRowsFiltered'); // wait to load monomers
      await awaitGrid(view.grid);
      //debugger;

      _package.logger.debug('Bio tests: substructureFilters/sync-helm, before changed event');
      await delay(f1.debounceTime * 2);
      _package.logger.debug('Bio tests: substructureFilters/sync-helm, after changed event');
      expect(df.filter.trueCount, fTrueCount);

      await f1.awaitRendered();
      expect(bf2.props.substructure, fSubStr);
    } finally {
      f1.detach();
      f2.detach();
      dlg.close();
    }
    await Promise.all([f1.awaitRendered(), f2.awaitRendered()]);
    await awaitGrid(view.grid);
    await delay(3000); //TODO: await for grid.onLookChanged
  });

  // two seq columns

  const twoColumnsCsv: string = `id,seq1,seq2,trueSeq1,trueSeq2
0,CGGCTACGGC,ATTGCATTCG,0,1,
1,CGGCTGCCGC,ATAGCATTCG,1,1,
2,CGGCTGCGCC,AATGCATACG,1,0,
3,CGGCTGCATT,TTTGCATTCG,1,1,
4,CGGCTGCATT,AAAGCATACG,1,0,
`;

  test('two-columns-fasta', async () => {
    const df = DG.DataFrame.fromCsv(twoColumnsCsv);
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);

    const fSeq1ColName: string = 'seq1';
    const fSeq1SubStr: string = 'CGGCTG';
    const fSeq1Trues: number[] = df.getCol('trueSeq1').toList();

    const fSeq2ColName: string = 'seq2';
    const fSeq2SubStr: string = 'GCATT';
    const fSeq2Trues: number[] = df.getCol('trueSeq2').toList();

    //const seq2Filter = new BioSubstructureFilter();
    const filterList: any[] = [
      {type: 'Bio:bioSubstructureFilter', columnName: fSeq1ColName},
      {type: 'Bio:bioSubstructureFilter', columnName: fSeq2ColName},
    ];
    const fg = (await df.plot.fromType(DG.VIEWER.FILTERS,
      {filters: filterList})) as DG.FilterGroup;
    view.dockManager.dock(fg, DG.DOCK_TYPE.LEFT);
    await delay(100);
    await awaitGrid(view.grid);

    const seq1Filter = fg.filters[0] as BioSubstructureFilter;
    const seq2Filter = fg.filters[1] as BioSubstructureFilter;
    expect(seq1Filter.column!.name, fSeq1ColName);
    expect(seq2Filter.column!.name, fSeq2ColName);

    const seq1Bf = seq1Filter.bioFilter as FastaBioFilter;
    const seq2Bf = seq2Filter.bioFilter as FastaBioFilter;

    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq1Bf.props = new BioFilterProps(fSeq1SubStr);
    }, 1000);
    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq2Bf.props = new BioFilterProps('');
    }, 1000, 'testEvent onRowsFiltered on seq1');
    expect(df.filter.trueCount, fSeq1Trues.filter((v) => v === 1).length);
    expect(df.filter.toBinaryString(), fSeq1Trues.map((v) => v.toString()).join(''));

    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq1Bf.props = new BioFilterProps('');
    }, 1000);
    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq2Bf.props = new BioFilterProps(fSeq2SubStr);
    }, 1000, 'testEvent onRowsFiltered on seq2');
    expect(df.filter.trueCount, fSeq2Trues.filter((v) => v === 1).length);
    expect(df.filter.toBinaryString(), fSeq2Trues.map((v) => v.toString()).join(''));

    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq1Bf.props = new BioFilterProps('');
    }, 1000);
    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq2Bf.props = new BioFilterProps('');
    }, 1000, 'testEvent onRowsFiltered on neither');
    expect(df.filter.trueCount, df.rowCount);

    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq1Bf.props = new BioFilterProps(fSeq1SubStr);
    }, 5000);
    await testEvent(df.onRowsFiltered, () => {}, () => {
      seq2Bf.props = new BioFilterProps(fSeq2SubStr);
    }, 5000, 'testEvent onRowsFiltered on both');
    const bothTrues: number[] = wu.count(0).take(df.rowCount)
      .map((rowI) => fSeq1Trues[rowI] * fSeq2Trues[rowI]).toArray();
    expect(df.filter.trueCount, bothTrues.filter((v) => v === 1).length);
    expect(df.filter.toBinaryString(), bothTrues.map((v) => v.toString()).join(''));

    await Promise.all([seq1Filter.awaitRendered(), seq2Filter.awaitRendered(), awaitGrid(view.grid)]);
  });

  // -- reset --

  test('reset-fasta', async () => {
    const df = await readDataframe('tests/filter_FASTA.csv');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);

    const fSeqColName: string = 'fasta';
    const fSubStr: string = 'MD';
    const fTrueCount: number = 3;

    const filterList = [{type: 'Bio:bioSubstructureFilter', columnName: fSeqColName}];
    const fg = (await df.plot.fromType(DG.VIEWER.FILTERS,
      {filters: filterList})) as DG.FilterGroup;
    view.dockManager.dock(fg, DG.DOCK_TYPE.LEFT);
    await delay(100);
    await awaitGrid(view.grid);

    const seqFilter = fg.filters[0] as BioSubstructureFilter;
    const seqBf = seqFilter.bioFilter as FastaBioFilter;
    await testEvent(df.onRowsFiltered, () => {}, () => {
      seqBf.props = new BioFilterProps(fSubStr);
    }, 1000, 'testEvent onRowsFiltered');
    expect(df.filter.trueCount, fTrueCount);
    expect(seqBf.props.substructure, fSubStr);
    expect(seqBf.substructureInput.value, fSubStr);

    const fgResetIconEl: HTMLElement = $(fg.root).find('i[name="icon-arrow-rotate-left"]')[0] as HTMLElement;
    fgResetIconEl.click();
    await delay(100);
    await awaitGrid(view.grid);
    expect(seqBf.props.substructure, '');
    expect(seqBf.substructureInput.value, '');
  });

  test('reopen', async () => {
    const df = await _package.files.readCsv('tests/filter_FASTA.csv');
    const view = grok.shell.addTableView(df);

    const filterList = [{type: 'Bio:bioSubstructureFilter', columnName: 'fasta'}];

    const fg1 = (await df.plot.fromType(DG.VIEWER.FILTERS,
      {filters: filterList})) as DG.FilterGroup;
    const fg1Dn = view.dockManager.dock(fg1, DG.DOCK_TYPE.LEFT);
    await delay(100);
    await awaitGrid(view.grid);
    fg1.close();
    await awaitGrid(view.grid);

    const fg2 = (await df.plot.fromType(DG.VIEWER.FILTERS,
      {filters: filterList})) as DG.FilterGroup;
    const fg2Dn = view.dockManager.dock(fg2, DG.DOCK_TYPE.LEFT);
    await delay(100);
    await awaitGrid(view.grid);
  });
});

async function createFilter(colName: string, df: DG.DataFrame): Promise<BioSubstructureFilter> {
  if (!df.columns.names().includes(colName)) {
    throw new Error(`The column '${colName}' not found. ` +
      `Available in data frame are ${JSON.stringify(df.columns.names())}`);
  }

  const filter = new BioSubstructureFilter();
  filter.attach(df);
  filter.applyState({columnName: colName});
  filter.column = df.col(colName);
  filter.columnName = colName;
  //filter.tableName = df.name;
  return filter;
};
