/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Bar Chart', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('barChart.visual', async () => {
    const barChartIcon = document.getElementsByClassName('svg-bar-chart')[0] as HTMLElement;
    barChartIcon.click();
    await awaitCheck(() => document.querySelector('.d4-bar-chart') !== null, 'bar chart not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Bar chart');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null, 'bar chart properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Bar chart viewer was not closed', 1000);
  });

  test('barChart.api', async () => {
    const barChart = v.barChart({
      split: 'race',
      value: 'height',
      valueAggrType: 'max',
    });
    await awaitCheck(() => document.querySelector('.d4-bar-chart') !== null, 'bar chart not found', 3000);

    if (barChart.props.splitColumnName != 'race')
      throw 'Split column has not been set'; 
    if (barChart.props.valueColumnName != 'height')
      throw 'Value column has not been set';
    if (barChart.props.valueAggrType != 'max')
      throw 'Aggregation has not been set'; 
    
    barChart.setOptions({
      title: 'Test Bar chart',
      relativeValues: true,
      onClick: 'Filter',
      barBorderLineMouseOverWidth: 10,
    });
    
    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Bar chart', 'title property has not been set', 2000);
    if (!barChart.props.relativeValues)
      throw 'relativeValues property has not been set'; 
    if (barChart.props.onClick != 'Filter')
      throw 'onClick property has not been set';
    if (barChart.props.barBorderLineMouseOverWidth !== 10)
      throw 'barBorderLineMouseOverWidth property has not been set';
    if (barChart.props.barSortOrder != 'desc')
      throw 'barSortOrder property (default) has not been set';
    if (barChart.props.barSortType != 'by value')
      throw 'barSortType property (default) has not been set';   
  });

  // Does not work through Test Manager
  test('barChart.serialization', async () => {
    await uploadProject('Test project with Bar Chart', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Bar Chart');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Bar chart');
    const barChart = findViewer('Bar chart', v);

    if (!barChart!.props.relativeValues)
      throw 'relativeValues property has not been deserialized'; 
    if (barChart!.props.onClick != 'Filter')
      throw 'onClick property has not been deserialized';
    if (barChart!.props.barBorderLineMouseOverWidth != 10)
      throw 'barBorderLineMouseOverWidth property has not been deserialized';
    if (barChart!.props.barSortOrder != 'desc')
      throw 'barSortOrder property (default) has not been deserialized';
    if (barChart!.props.barSortType != 'by value')
      throw 'barSortType property (default) has not been deserialized';
    if (barChart!.props.splitColumnName != 'race')
      throw 'Split column has not been deserialized'; 
    if (barChart!.props.valueColumnName != 'height')
      throw 'Value column has not been deserialized';
    if (barChart!.props.valueAggrType != 'max')
      throw 'Aggregation has not been deserialized';
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Bar Chart').first());
  }); 
});
