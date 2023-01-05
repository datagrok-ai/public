/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Scatter Plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('scatterPlot.visual', async () => {
    const scatterPlotIcon = document.getElementsByClassName('svg-scatter-plot')[0] as HTMLElement;
    scatterPlotIcon.click(); 
    await awaitCheck(() => document.querySelector('.d4-scatter-plot') !== null, 'scatter plot not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null,
      'scatter plot properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Scatter plot viewer was not closed', 1000);

    /* finish when it becomes clear with menu events
    let markerMenu:HTMLElement;
    let contextMenuItem;
    for (let i=0; i<document.getElementsByClassName('d4-menu-item-label').length; i++) {
        contextMenuItem = document.getElementsByClassName('d4-menu-item-label')[i] as HTMLElement;
        if (contextMenuItem.innerText == 'Markers'){
            markerMenu = contextMenuItem;
            break;
        }
    } */
  });

  test('scatterPlot.api', async () => {
    const scatterPlot = v.scatterPlot({
      x: 'weight',
      y: 'height',
      size: 'age',
      color: 'sex',
    });
    await awaitCheck(() => document.querySelector('.d4-scatter-plot') !== null, 'scatter plot not found', 3000);

    if (scatterPlot.props.xColumnName != 'weight')
      throw 'X column has not been set'; 
    if (scatterPlot.props.yColumnName != 'height')
      throw 'Y column has not been set';
    if (scatterPlot.props.sizeColumnName != 'age')
      throw 'Size column has not been set'; 
    if (scatterPlot.props.colorColumnName != 'sex')
      throw 'Name column has not been set';

    scatterPlot.setOptions({
      title: 'Test Scatter Plot',
      showRegressionLine: true,
      markerType: 'dot',
      showVerticalGridLines: true,
      markerBorderWidth: 10, 
    });

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Scatter Plot', 'title property has not been set', 2000);
    if (!scatterPlot.props.showRegressionLine)
      throw 'showRegressionLine property has not been set'; 
    if (scatterPlot.props.markerType != 'dot')
      throw 'markerType property has not been set';
    if (!scatterPlot.props.showVerticalGridLines)
      throw 'showVerticalGridLines property has not been set';
    if (scatterPlot.props.markerBorderWidth != 10)
      throw 'markerBorderWidth property has not been set';
  });

  // Does not work through Test Manager
  test('scatterPlot.serialization', async () => {
    await uploadProject('Test project with Scatter plot', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Scatter plot');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    const scatterPlot = findViewer('Scatter plot', v);

    if (scatterPlot!.props.xColumnName != 'weight')
      throw 'X column has not been deserialized'; 
    if (scatterPlot!.props.yColumnName != 'height')
      throw 'Y column has not been deserialized';
    if (scatterPlot!.props.sizeColumnName != 'age')
      throw 'Size column has not been deserialized'; 
    if (scatterPlot!.props.colorColumnName != 'sex')
      throw 'Name column has not been deserialized';
    if (!scatterPlot!.props.showRegressionLine)
      throw 'showRegressionLine property has not been deserialized'; 
    if (scatterPlot!.props.markerType != 'dot')
      throw 'markerType property has not been deserialized';
    if (!scatterPlot!.props.showVerticalGridLines)
      throw 'showVerticalGridLines property has not been deserialized';
    if (scatterPlot!.props.markerBorderWidth != 10)
      throw 'markerBorderWidth property has not been deserialized';         
  });
  
  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Scatter plot').first());
  }); 
});
