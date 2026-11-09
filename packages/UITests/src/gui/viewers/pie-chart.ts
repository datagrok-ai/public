/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, settled, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Pie Chart', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
  });

  test('pieChart.visual', async () => {
    v = grok.shell.addTableView(demog);
    const pieChartIcon = document.getElementsByClassName('svg-pie-chart')[0] as HTMLElement;
    pieChartIcon.click(); 
    await waitUntil(() => document.querySelector('.d4-pie-chart') !== null, 'pie chart not found');
    isViewerPresent(Array.from(v.viewers), 'Pie chart');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await waitUntil(() => document.querySelector('.grok-prop-panel') !== null, 'pie chart properties not found');
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await waitUntil(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found');
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await waitUntil(() => Array.from(v.viewers).length === 1, 'Pie chart viewer was not closed');
  });

  test('pieChart.api', async () => {
    v = grok.shell.addTableView(demog);
    const pieChart = v.pieChart({
      category: 'race',
    });
    await waitUntil(() => document.querySelector('.d4-pie-chart') !== null, 'pie chart not found');

    if (pieChart.props.categoryColumnName != 'race')
      throw 'Category column has not been set'; 
    
    pieChart.setOptions({
      title: 'Test Pie chart',
      includeNulls: false,
      segmentLengthColumnName: 'age',
      onClick: 'Filter',
      startAngle: 30, 
    });

    await waitUntil(() => (Array.from(document.querySelectorAll(
      '#elementContent > div.d4-layout-top > div > textarea')) as HTMLTextAreaElement[])
      .some((e) => e.value === 'Test Pie chart'), 'title property has not been set');
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
    v = grok.shell.addTableView(demog);
    await settled(v.pieChart({
      category: 'race',
      title: 'Test Pie chart',
      includeNulls: false,
      segmentLengthColumnName: 'age',
      onClick: 'Filter',
      startAngle: 30,
    }));
    const project = await uploadProject('Test project with Pie Chart', demog.getTableInfo(), v, demog);
    try {
      grok.shell.closeAll();
      await project.open();
      v = grok.shell.tableView('demog 1000');
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
    } finally {
      await grok.dapi.projects.delete(project);
    }
  }); 
}, { owner: 'dkovalyov@datagrok.ai' });
