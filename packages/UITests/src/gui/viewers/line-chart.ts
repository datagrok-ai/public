/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {isViewerPresent, findViewer, showToolbox} from '../gui-utils';


category('Viewers: Line Chart', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    showToolbox();
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('lineChart.visual', async () => {
    grok.shell.windows.showToolbox = true;
    await delay(100);
    const lineChartIcon = document.getElementsByClassName('svg-line-chart')[0] as HTMLElement;
    lineChartIcon.click();
    await awaitCheck(() => document.querySelector('.d4-line-chart') !== null, 'line chart not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Line chart');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null,
      'line chart properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Line chart viewer was not closed', 1000);
  });

  test('lineChart.api', async () => {
    const lineChart = v.lineChart({
      xColumnName: 'age',
      yColumnNames: ['height', 'weight'],
      yAggrTypes: ['avg', 'min'],
      splitColumnName: 'race',
    });
    await awaitCheck(() => document.querySelector('.d4-line-chart') !== null, 'line chart not found', 3000);

    if (lineChart.props.splitColumnName != 'race')
      throw 'Split column has not been set';
    if (lineChart.props.yColumnNames[0] != 'height')
      throw 'Value column has not been set';
    if (lineChart.props.yAggrTypes[0] != 'avg')
      throw 'Aggregation has not been set';

    lineChart.setOptions({
      title: 'Test Line chart',
    });

    await awaitCheck(() => (document.querySelector('.d4-line-chart')!.parentElement!.parentElement!
      .querySelector('.panel-titlebar') as HTMLElement)!.innerText === 'Test Line chart', 'title property has not been set', 2000);
  });

  // Verifies that Line Chart settings survive layout (view state) serialization —
  // the same serialization the platform uses to persist views in projects.
  test('lineChart.serialization', async () => {
    const source = grok.shell.addTableView(grok.data.demo.demog(1000));
    const lineChart = source.lineChart({
      xColumnName: 'age',
      yColumnNames: ['height', 'weight'],
      yAggrTypes: ['avg', 'min'],
      splitColumnName: 'race',
    });
    await awaitCheck(() => document.querySelector('.d4-line-chart') !== null, 'line chart not found', 3000);
    lineChart.setOptions({lineWidth: 3, showMarkers: 'Always', multiAxis: true, invertXAxis: true});
    const layout = source.saveLayout();

    const target = grok.shell.addTableView(grok.data.demo.demog(1000));
    target.loadLayout(layout);
    await awaitCheck(() => findViewer('Line chart', target) !== undefined,
      'Line chart not found after layout deserialization', 3000);
    isViewerPresent(Array.from(target.viewers), 'Line chart');
    const restored = findViewer('Line chart', target);

    if (restored!.props.splitColumnName != 'race')
      throw 'Split column has not been deserialized';
    if (restored!.props.xColumnName != 'age')
      throw 'X column has not been deserialized';
    if (restored!.props.yColumnNames[0] != 'height')
      throw 'Value column has not been deserialized';
    if (restored!.props.yColumnNames[1] != 'weight')
      throw 'Second value column has not been deserialized';
    if (restored!.props.yAggrTypes[0] != 'avg')
      throw 'Aggregation has not been deserialized';
    if (restored!.props.lineWidth != 3)
      throw 'lineWidth property has not been deserialized';
    if (restored!.props.showMarkers != 'Always')
      throw 'showMarkers property has not been deserialized';
    if (!restored!.props.multiAxis)
      throw 'multiAxis property has not been deserialized';
    if (!restored!.props.invertXAxis)
      throw 'invertXAxis property has not been deserialized';
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, { owner: 'dkovalyov@datagrok.ai' });
