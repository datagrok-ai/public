/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Pie Chart', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('pieChart.visual', async () => {
    const pieChartIcon = document.getElementsByClassName('svg-pie-chart')[0] as HTMLElement;
    pieChartIcon.click(); 
    await awaitCheck(() => document.querySelector('.d4-pie-chart') !== null, 'pie chart not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Pie chart');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null, 'pie chart properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Pie chart viewer was not closed', 1000);
  });

  test('pieChart.api', async () => {
    const pieChart = v.pieChart({
      category: 'race',
    });
    await awaitCheck(() => document.querySelector('.d4-pie-chart') !== null, 'pie chart not found', 3000);

    if (pieChart.props.categoryColumnName != 'race')
      throw 'Category column has not been set'; 
    
    pieChart.setOptions({
      title: 'Test Pie chart',
      includeNulls: false,
      segmentLengthColumnName: 'age',
      onClick: 'Filter',
      startAngle: 30, 
    });

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Pie chart', 'title property has not been set', 2000);
    if (pieChart.props.includeNulls)
      throw 'includeNulls property has not been set to false'; 
    if (pieChart.props.segmentLengthColumnName != 'age')
      throw 'segmentLengthColumnName property has not been set';
    if (pieChart.props.onClick != 'Filter')
      throw 'onClick property has not been set';
    if (pieChart.props.startAngle != 30)
      throw 'startAngle property (default) has not been set';
  });  

  // Does not work through Test Manager
  test('pieChart.serialization', async () => {
    await uploadProject('Test project with Pie Chart', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Pie Chart');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Pie chart');
    const pieChart = findViewer('Pie chart', v);

    if (pieChart!.props.includeNulls)
      throw 'includeNulls property has not been deserialized'; 
    if (pieChart!.props.segmentLengthColumnName != 'age')
      throw 'segmentLengthColumnName property has not been deserialized';
    if (pieChart!.props.onClick != 'Filter')
      throw 'onClick property has not been deserialized';
    if (pieChart!.props.startAngle != 30)
      throw 'startAngle property (default) has not been deserialized';
    if (pieChart!.props.categoryColumnName != 'race')
      throw 'Category column has not been deserialized'; 
  }); 
  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Pie Chart').first());
  }); 
});
