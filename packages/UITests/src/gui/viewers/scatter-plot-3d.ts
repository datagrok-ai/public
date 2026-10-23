/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, settled, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: 3D Scatter Plot', () => {  
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
  });

  test('scatterPlot3D.visual', async () => {
    v = grok.shell.addTableView(demog);
    const ScatterPlot3DIcon = document.getElementsByClassName('svg-3d-scatter-plot')[0] as HTMLElement;
    ScatterPlot3DIcon.click();
    await waitUntil(() => document.querySelector('.d4-3d-scatter-plot') !== null, '3D scatter plot not found');
    isViewerPresent(Array.from(v.viewers), '3d scatter plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await waitUntil(() => document.querySelector('.grok-prop-panel') !== null,
      '3D scatter plot properties not found');
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await waitUntil(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found');
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await waitUntil(() => Array.from(v.viewers).length === 1, '3D scatter plot viewer was not closed');
  });

  test('scatterPlot3D.api', async () => {
    v = grok.shell.addTableView(demog);
    const scatterPlot3D = v.addViewer(DG.VIEWER.SCATTER_PLOT_3D, {
      'xColumnName': 'height',
      'yColumnName': 'weight',
      'zColumnName': 'age',
      'sizeColumnName': 'race',
      'colorColumnName': 'sex',
    });
    await waitUntil(() => document.querySelector('.d4-3d-scatter-plot') !== null, '3D scatter plot not found');

    if (scatterPlot3D.props.xColumnName != 'height')
      throw 'X column has not been set'; 
    if (scatterPlot3D.props.yColumnName != 'weight')
      throw 'Y column has not been set';
    if (scatterPlot3D.props.zColumnName != 'age')
      throw 'Z column has not been set'; 
    if (scatterPlot3D.props.sizeColumnName != 'race')
      throw 'Size column has not been set';
    if (scatterPlot3D.props.colorColumnName != 'sex')
      throw 'Color column has not been set'; 
    
    scatterPlot3D.setOptions({
      title: 'Test 3D Scatter Plot',
      dynamicCameraMovement: true,
      markerType: 'sphere',
      markerOpacity: 17,
    });

    await waitUntil(() => scatterPlot3D.props.title === 'Test 3D Scatter Plot',
      'title property has not been set');  
    if (!scatterPlot3D.props.dynamicCameraMovement)
      throw 'dynamicCameraMovement property has not been set'; 
    if (scatterPlot3D.props.markerType != 'sphere')
      throw 'markerType property has not been set';
    if (scatterPlot3D.props.markerOpacity != 17)
      throw 'markerOpacity property has not been set';
  });

  // Does not work through Test Manager
  test('scatterPlot3D.serialization', async () => {
    v = grok.shell.addTableView(demog);
    await settled(v.addViewer(DG.VIEWER.SCATTER_PLOT_3D, {
      xColumnName: 'height',
      yColumnName: 'weight',
      zColumnName: 'age',
      sizeColumnName: 'race',
      colorColumnName: 'sex',
      title: 'Test 3D Scatter Plot',
      dynamicCameraMovement: true,
      markerType: 'sphere',
      markerOpacity: 17,
    }));
    const project = await uploadProject('Test project with 3 D Scatter Plot', demog.getTableInfo(), v, demog);
    try {
      grok.shell.closeAll();
      await project.open();
      v = grok.shell.tableView('demog 1000');
      isViewerPresent(Array.from(v.viewers), '3d scatter plot');
      const scatterPlot3D = findViewer('3d scatter plot', v);

      if (!scatterPlot3D!.props.dynamicCameraMovement)
        throw 'dynamicCameraMovement property has not been deserialized'; 
      if (scatterPlot3D!.props.markerType != 'sphere')
        throw 'markerType property has not been deserialized';
      if (scatterPlot3D!.props.markerOpacity != 17)
        throw 'markerOpacity property has not been deserialized';
      if (scatterPlot3D!.props.xColumnName != 'height')
        throw 'X column has not been deserialized'; 
      if (scatterPlot3D!.props.yColumnName != 'weight')
        throw 'Y column has not been deserialized';
      if (scatterPlot3D!.props.zColumnName != 'age')
        throw 'Z column has not been deserialized'; 
      if (scatterPlot3D!.props.sizeColumnName != 'race')
        throw 'Size column has not been deserialized';
      if (scatterPlot3D!.props.colorColumnName != 'sex')
        throw 'Color column has not been deserialized'; 
    } finally {
      await grok.dapi.projects.delete(project);
    }
  }); 
}, { owner: 'dkovalyov@datagrok.ai' });
