/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';

  
category('Viewers: Box Plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('boxPlot.visual', async () => {
    const boxPlotIcon = document.getElementsByClassName('svg-box-plot')[0] as HTMLElement;
    boxPlotIcon.click();
    await awaitCheck(() => document.querySelector('.d4-box-plot') !== null, 'box plot not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Box plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null, 'box plot properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Box plot viewer was not closed', 1000);
  }); 

  test('boxPlot.api', async () => {
    const boxPlot = v.addViewer(DG.VIEWER.BOX_PLOT, {
      'valueColumnName': 'age',
      'categoryColumnName': 'race',
      'binColorColumnName': 'height',
      'markerColorColumnName': 'sex',
    });
    await awaitCheck(() => document.querySelector('.d4-box-plot') !== null, 'box plot not found', 3000);

    if (boxPlot.props.valueColumnName != 'age')
      throw 'Value column has not been set'; 
    if (boxPlot.props.categoryColumnName != 'race')
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

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Box Plot', 'title property has not been set', 2000);
    if (boxPlot.props.markerSize != 10)
      throw 'Marker size property has not been set to 10'; 
    if (!boxPlot.props.showStatistics)
      throw 'Show Statistics property has not been set to TRUE'; 
    if (boxPlot.props.axisType != 'logarithmic')
      throw 'Axis Type property has not been set to logarithmic';
  });

  // Does not work through Test Manager
  test('boxPlot.serialization', async () => {
    await uploadProject('Test project with Box Plot', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Box Plot');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Box plot');
    const boxPlot = findViewer('Box plot', v);

    if (boxPlot!.props.valueColumnName != 'age')
      throw 'Value column has not been deserialized'; 
    if (boxPlot!.props.categoryColumnName != 'race')
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
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Box Plot').first());
  });
});
