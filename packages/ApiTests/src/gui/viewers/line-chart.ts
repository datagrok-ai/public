/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Line Chart', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('lineChart.visual', async () => {
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
    });
    await awaitCheck(() => document.querySelector('.d4-line-chart') !== null, 'line chart not found', 3000);

    if (lineChart.props.splitColumnName != 'race')
      throw 'Split column has not been set'; 
    if (lineChart.props.valueColumnName != 'height')
      throw 'Value column has not been set';
    if (lineChart.props.valueAggrType != 'max')
      throw 'Aggregation has not been set'; 

    lineChart.setOptions({
      title: 'Test Line chart',
      relativeValues: true,
      onClick: 'Filter',
      barBorderLineMouseOverWidth: 10,
    });

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Line chart', 'title property has not been set', 2000);
    if (!lineChart.props.relativeValues)
      throw 'relativeValues property has not been set'; 
    if (lineChart.props.onClick != 'Filter')
      throw 'onClick property has not been set';
    if (lineChart.props.barBorderLineMouseOverWidth != 10)
      throw 'barBorderLineMouseOverWidth property has not been set';
    if (lineChart.props.barSortOrder != 'desc')
      throw 'barSortOrder property (default) has not been set';
    if (lineChart.props.barSortType != 'by value')
      throw 'barSortType property (default) has not been set';   
  }, {skipReason: 'need fix'});

  // Does not work through Test Manager
  test('lineChart.serialization', async () => {
    await uploadProject('Test project with Line Chart', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Line Chart');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Line chart');
    const lineChart = findViewer('Line chart', v);

    if (lineChart!.props.splitColumnName != 'race')
      throw 'Split column has not been set'; 
    if (lineChart!.props.valueColumnName != 'height')
      throw 'Value column has not been set';
    if (lineChart!.props.valueAggrType != 'max')
      throw 'Aggregation has not been set';
    if (!lineChart!.props.relativeValues)
      throw 'relativeValues property has not been set'; 
    if (lineChart!.props.onClick != 'Filter')
      throw 'onClick property has not been set';
    if (lineChart!.props.barBorderLineMouseOverWidth != 10)
      throw 'barBorderLineMouseOverWidth property has not been set';
    if (lineChart!.props.barSortOrder != 'desc')
      throw 'barSortOrder property (default) has not been set';
    if (lineChart!.props.barSortType != 'by value')
      throw 'barSortType property (default) has not been set';   
  }, {skipReason: 'need previous test'});

  after(async () => {
    grok.shell.closeAll();
    // await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Line Chart').first());
  }); 
});
