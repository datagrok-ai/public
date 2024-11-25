/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Histogram', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
  });

  test('histogram.visual', async () => {
    v = grok.shell.addTableView(demog);
    const histogramChartIcon = document.getElementsByClassName('svg-histogram')[0] as HTMLElement;
    histogramChartIcon.click();
    await awaitCheck(() => document.querySelector('.d4-histogram') !== null, 'histogram not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Histogram');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null, 'histogram properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Histogram viewer was not closed', 1000);
  });

  test('histogram.api', async () => {
    v = grok.shell.addTableView(demog);
    const histogram = v.histogram({
      value: 'weight',
      binWidthRatio: 2,
    });
    await awaitCheck(() => document.querySelector('.d4-histogram') !== null, 'histogram not found', 3000);

    if (histogram.props.valueColumnName != 'weight')
      throw 'Value column has not been set';
    if (histogram.props.binWidthRatio != 2)
      throw 'binWidthRatio has not been set to 2';
    const filterCheckbox = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0]
      .getElementsByClassName('d4-filter-out-missing-values')[0] as HTMLInputElement;
    if (!filterCheckbox)
      throw 'filterCheckbox element not found';

    histogram.setOptions({
      title: 'Test Histogram',
      colorColumnName: 'age',
      showRangeInputs: true,
      bins: 30,
    });

    await awaitCheck(() => (Array.from(document.querySelectorAll('#elementContent > div.d4-layout-top > div > textarea')) as HTMLSelectElement[]).some((e)=>e.value === 'Test Histogram')
      , 'title property has not been set', 2000);
    if (histogram.props.bins != 30)
      throw 'bins property has not been set to 30';
    if (histogram.props.colorColumnName != 'age')
      throw 'colorColumnName property has not been set';
    const rangeInputMin = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0]
      .getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-min')[0] as HTMLInputElement;
    const rangeInputMax = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0]
      .getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-max')[0] as HTMLInputElement;
    if (!rangeInputMin)
      throw 'min range input from showRangeInputs property not found';
    if (!rangeInputMax)
      throw 'max range input from showRangeInputs property not found';
  });

  // Does not work through Test Manager
  test('histogram.serialization', async () => {
    v = grok.shell.addTableView(demog);
    await uploadProject('Test project with Histogram', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Histogram');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Histogram');
    const histogram = findViewer('Histogram', v);

    if (histogram!.props.valueColumnName != 'weight')
      throw 'Value column has not been deserialized';
    if (histogram!.props.binWidthRatio != 2)
      throw 'binWidthRatio has not been deserialized';
    if (histogram!.props.bins != 30)
      throw 'bins property has not been deserialized';
    if (histogram!.props.colorColumnName != 'age')
      throw 'colorColumnName property has not been deserialized';
    if (!histogram!.props.showRangeInputs)
      throw 'showRangeInputs property has not been deserialized';

    const filterCheckbox = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0]
      .getElementsByClassName('d4-filter-out-missing-values')[0] as HTMLInputElement;
    if (!filterCheckbox)
      throw 'filterCheckbox element not found';
    const rangeInputMin = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0]
      .getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-min')[0] as HTMLInputElement;
    const rangeInputMax = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0]
      .getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-max')[0] as HTMLInputElement;
    if (!rangeInputMin)
      throw 'min range input from showRangeInputs property not found';
    if (!rangeInputMax)
      throw 'max range input from showRangeInputs property not found';
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Histogram').first());
  }, {skipReason: 'GROK-12698'});

  test('histogram.spline.range.overlapped', async () => {
    const df: DG.DataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [0, 6, 8])]);
    df.getCol('value').setTag('.charts', JSON.stringify([{
      title: 'Test spline', type: 'spline',
      color: '#7570B3', width: 1, ['normalize-y']: true, visible: true,
      x: [2, 7, 10],
      y: [1, 2, 1],
    }]));
    const tv = grok.shell.addTableView(df);

    try {
      const viewer = (await df.plot.fromType(DG.VIEWER.HISTOGRAM, {})) as DG.Viewer;
      tv.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Histogram test', 0.3);

      await awaitCheck(() => document.querySelector('.d4-histogram') !== null, 'histogram not found', 3000);
      isViewerPresent(Array.from(tv.viewers), 'Histogram');
    } finally {
      tv.close();
      grok.shell.closeTable(df);
    }
  });

  test('histogram.spline.range.non-overlapped', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', [1, 2, 3]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [1000, 1006, 1008])]);
    df.getCol('value').setTag('.charts', JSON.stringify([{
      title: 'Test spline', type: 'spline',
      color: '#7570B3', width: 1, ['normalize-y']: true, visible: true,
      x: [1, 3, 12],
      y: [1, 3, 2],
    }]));
    const tv = grok.shell.addTableView(df);

    try {
      const viewer = (await df.plot.fromType(DG.VIEWER.HISTOGRAM, {value: 'value'})) as DG.Viewer;
      tv.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT, null, 'Histogram test', 0.3);

      await awaitCheck(() => document.querySelector('.d4-histogram') !== null, 'histogram not found', 3000);
      isViewerPresent(Array.from(tv.viewers), 'Histogram');
    } finally {
      tv.close();
      grok.shell.closeTable(df);
    }
  });

  after(async () => {
    grok.shell.closeAll();
  });
});
