// import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { category, test, before, after, settled, waitUntil } from '@datagrok-libraries/test/src/test';
import $ from 'cash-dom';

category('Viewers: Filters', () => {
  const csv1: string = `id1,id2,id3
id1_0001,id2_001,id3_1
id1_0002,id2_002,id3_2
id1_0003,id2_003,id3_3`;

  test('twoCategorical', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv1);
    const view: DG.TableView = grok.shell.addTableView(df);

    const filterList: { [p: string]: string }[] = [
      { column: 'id1', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id1 label' },
      { column: 'id2', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id2 label' },
    ];
    const filtersViewer = view.filters({ filters: filterList }) as DG.Viewer;
    view.dockManager.dock(filtersViewer, DG.DOCK_TYPE.LEFT, null, 'Filters', 0.4);

    await waitUntil(() => $(filtersViewer.root)
      .find('div.d4-flex-row.d4-filter-header label.d4-filter-column-name').get()
      .map((lbl) => lbl.innerText).join() === ['id1', 'id2'].join(), 'cannot find all filters');
  });

  test('customBetweenCategorical', async () => {
    let _package = DG.Func.find({ package: 'UITests', name: 'test' })[0]?.package;
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv1);
    const view: DG.TableView = grok.shell.addTableView(df);

    const filterList: { [p: string]: string }[] = [
      { column: 'id1', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id1 label' },
      { column: 'id3', type: `${_package.name}:TestCustomFilter`, label: 'custom label' },
      { column: 'id2', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id2 label' },
    ];
    const filtersViewer = view.filters({ filters: filterList }) as DG.Viewer;

    await waitUntil(() => $(filtersViewer.root)
      .find('div.d4-flex-row.d4-filter-header label.d4-filter-column-name').get()
      .map((lbl) => lbl.innerText).join() === ['id1', 'id3', 'id2'].join(), 'cannot find all filters');
  });
}, { owner: 'drizhinashvili@datagrok.ai' });


category('Viewers: Filters: Collaborative filtering', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;
  let fg: DG.FilterGroup;
  let _package: DG.Package;
  const STRUCTURE = 'Structure';
  const BENZENE = 'c1ccccc1';

  before(async () => {
    _package = DG.Func.find({ package: 'UITests', name: 'test' })[0]?.package;
    await DG.Func.find({package: 'Chem', name: 'scaffoldTreeFilter'})[0].package.load();
    df = await _package.files.readCsv('SPGI_v2_100_full.csv');
    tv = grok.shell.addTableView(df);
    await waitUntil(() => document.querySelector('canvas') !== null, 'cannot load table');
    fg = tv.getFiltersGroup();
  });

  test('Scatter plot', async () => {
    await addScatterPlot();
    await addScaffoldFilter(47);
    await addSubscructureFilter(17);
    await addCategoricalFilter(9);
    await addHistogtramFilter(6);
    tv.resetLayout();
  });

  // FILTERS

  async function addCategoricalFilter(n: number = 36) {
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.CATEGORICAL,
      column: 'Stereo Category',
      selected: ['R_ONE'],
    });
    await waitUntil(() => df.filter.trueCount === n,
      () => `Categorical filter: expected ${n} rows, got ${df.filter.trueCount}`);
  }

  async function addHistogtramFilter(n: number = 80) {
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.HISTOGRAM,
      column: 'Average Mass',
      min: 300,
      max: 460,
    });
    await waitUntil(() => df.filter.trueCount === n,
      () => `Histogram filter: expected ${n} rows, got ${df.filter.trueCount}`);
  }

  async function addSubscructureFilter(n: number = 32) {
    fg.updateOrAdd({
      type: DG.FILTER_TYPE.SUBSTRUCTURE,
      column: STRUCTURE,
      columnName: STRUCTURE,
      molBlock: BENZENE,
    });
    await waitUntil(() => df.filter.trueCount === n,
      () => `Subscructure filter: expected ${n} rows, got ${df.filter.trueCount}`);
  }

  async function addScaffoldFilter(n: number = 92) {
    fg.updateOrAdd({
      type: 'Chem:scaffoldTreeFilter',
      column: STRUCTURE,
      columnName: STRUCTURE,
      savedTree,
    });
    await waitUntil(() => df.filter.trueCount === n,
      () => `Scaffold Tree filter: expected ${n} rows, got ${df.filter.trueCount}`);
  }

  // VIEWERS

  async function addScatterPlot(n: number = 49) {
    const v = DG.Viewer.scatterPlot(df, {
      xColumnName: 'Chemical Space X',
      yColumnName: 'Chemical Space Y',
      axesFollowFilter: false,
      zoomAndFilter: 'filter by zoom',
    });
    tv.addViewer(v);
    await settled(v);
    v.zoom(-5, 4, 3, 11);
    await waitUntil(() => df.filter.trueCount === n,
      () => `Scatter plot: expected ${n} rows, got ${df.filter.trueCount}`);
  }
}, { clear: false });

const savedTree = '[{"scaffold":"\\n     RDKit          2D\\n\\n  3  2  0  0  0  0  0  0  0  0999 V2000\\n    ' +
  '0.2933    2.5312    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\\n   -1.0463    1.7578    0.0000 C   ' +
  '0  0  0  0  0  0  0  0  0  0  0  0\\n    0.2933    4.0781    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  ' +
  '0\\n  2  1  1  0\\n  1  3  1  0\\nM  END\\n","checked":true,"isNot":false,"expanded":true,"child_nodes":[]},"OR"]';
