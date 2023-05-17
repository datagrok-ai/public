import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  category,
  test,
  awaitCheck,
  expectArray,
} from '@datagrok-libraries/utils/src/test';
import {TableView} from 'datagrok-api/dg';
import {_package} from '../package-test';

import $ from 'cash-dom';

category('Viewers: Filters', () => {
  const csv1: string = `id1,id2,id3
id1_0001,id2_001,id3_1
id1_0002,id2_002,id3_2
id1_0003,id2_003,id3_3`;

  test('twoCategorical', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv1);
    const view: TableView = grok.shell.addTableView(df);

    const filterList: { [p: string]: string }[] = [
      {column: 'id1', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id1 label'},
      {column: 'id2', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id2 label'},
    ];
    const filtersViewer = view.filters({filters: filterList}) as DG.Viewer;
    view.dockManager.dock(filtersViewer, DG.DOCK_TYPE.LEFT, null, 'Filters', 0.4);

    await awaitCheck(() => {
      const fltColNameList = $(filtersViewer.root)
        .find('div.d4-flex-row.d4-filter-header label.d4-filter-column-name').get()
        .map((lbl) => lbl.innerText);

      expectArray(fltColNameList, ['id1', 'id2']);
      return true;
    }, 'cannot find all filters', 3000);
  });

  test('customBetweenCategorical', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv1);
    const view: TableView = grok.shell.addTableView(df);

    const filterList: { [p: string]: string }[] = [
      {column: 'id1', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id1 label'},
      {column: 'id3', type: `${_package.name}:TestCustomFilter`, label: 'custom label'},
      {column: 'id2', type: DG.FILTER_TYPE.CATEGORICAL, label: 'id2 label'},
    ];
    const filtersViewer = view.filters({filters: filterList}) as DG.Viewer;

    await awaitCheck(() => {
      const fltColNameList = $(filtersViewer.root)
        .find('div.d4-flex-row.d4-filter-header label.d4-filter-column-name').get()
        .map((lbl) => lbl.innerText);

      expectArray(fltColNameList, ['id1', 'id3', 'id2']);
      return true;
    }, 'cannot find all filters', 3000);
  });
});
