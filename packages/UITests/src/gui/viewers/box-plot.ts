/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, settled, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';

  
category('Viewers: Box Plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
  });

  test('boxPlot.visual', async () => {
    v = grok.shell.addTableView(demog);
    const boxPlotIcon = document.getElementsByClassName('svg-box-plot')[0] as HTMLElement;
    boxPlotIcon.click();
    await waitUntil(() => document.querySelector('.d4-box-plot') !== null, 'box plot not found');
    isViewerPresent(Array.from(v.viewers), 'Box plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await waitUntil(() => document.querySelector('.grok-prop-panel') !== null, 'box plot properties not found');
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await waitUntil(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found');
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await waitUntil(() => Array.from(v.viewers).length === 1, 'Box plot viewer was not closed');
  }); 

  test('boxPlot.api', async () => {
    v = grok.shell.addTableView(demog);
    const boxPlot = v.addViewer(DG.VIEWER.BOX_PLOT, {
      'valueColumnName': 'age',
      'categoryColumnNames': ['race'],
      'binColorColumnName': 'height',
      'markerColorColumnName': 'sex',
    });
    await waitUntil(() => document.querySelector('.d4-box-plot') !== null, 'box plot not found');

    if (boxPlot.props.valueColumnName != 'age')
      throw 'Value column has not been set'; 
    if (boxPlot.props.categoryColumnNames[0] != 'race')
      throw 'Category column has not been set';     
    if (boxPlot.props.binColorColumnName != 'height')
      throw 'Bin color column has not been set';   
    if (boxPlot.props.markerColorColumnName != 'sex')
      throw 'Marker color column has not been set';    

    boxPlot.setOptions({
      title: 'Test Box Plot',
      markerSize: 10,
      showStatistics: true,
      axisType: 'logarithmic',
    });

    await waitUntil(() => (Array.from(document.querySelectorAll(
      '#elementContent > div.d4-layout-top > div > textarea')) as HTMLTextAreaElement[])
      .some((e) => e.value === 'Test Box Plot'), 'title property has not been set');
    if (boxPlot.props.markerSize != 10)
      throw 'Marker size property has not been set to 10'; 
    if (!boxPlot.props.showStatistics)
      throw 'Show Statistics property has not been set to TRUE'; 
    if (boxPlot.props.axisType != 'logarithmic')
      throw 'Axis Type property has not been set to logarithmic';
  });

  // Does not work through Test Manager
  test('boxPlot.serialization', async () => {
    v = grok.shell.addTableView(demog);
    await settled(v.addViewer(DG.VIEWER.BOX_PLOT, {
      valueColumnName: 'age',
      categoryColumnNames: ['race'],
      binColorColumnName: 'height',
      markerColorColumnName: 'sex',
      title: 'Test Box Plot',
      markerSize: 10,
      showStatistics: true,
      axisType: 'logarithmic',
    }));
    const project = await uploadProject('Test project with Box Plot', demog.getTableInfo(), v, demog);
    try {
      grok.shell.closeAll();
      await project.open();
      v = grok.shell.tableView('demog 1000');
      isViewerPresent(Array.from(v.viewers), 'Box plot');
      const boxPlot = findViewer('Box plot', v);

      if (boxPlot!.props.valueColumnName != 'age')
        throw 'Value column has not been deserialized'; 
      if (boxPlot!.props.categoryColumnNames[0] != 'race')
        throw 'Category column has not been deserialized';     
      if (boxPlot!.props.binColorColumnName != 'height')
        throw 'Bin color column has not been deserialized';   
      if (boxPlot!.props.markerColorColumnName != 'sex')
        throw 'Marker color column has not been deserialized';   
      if (boxPlot!.props.markerSize != 10)
        throw 'Marker size property has not been deserialized'; 
      if (!boxPlot!.props.showStatistics)
        throw 'Show Statistics property hhas not been deserialized'; 
      if (boxPlot!.props.axisType != 'logarithmic')
        throw 'Axis Type property has not been deserialized';
    } finally {
      await grok.dapi.projects.delete(project);
    }
  });
}, { owner: 'dkovalyov@datagrok.ai' });
