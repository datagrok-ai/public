/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Matrix Plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('matrixPlot.visual', async () => {
    const matrixPlotIcon = document.getElementsByClassName('svg-matrix-plot')[0] as HTMLElement;
    matrixPlotIcon.click(); 
    await awaitCheck(() => document.querySelector('.d4-matrix-plot') !== null, 'matrix plot not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Matrix plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null,
      'matrix plot properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Matrix plot viewer was not closed', 1000);
  });

  test('matrixPlot.api', async () => {
    const matrixPlot = v.matrixPlot({
      xColumnNames: ['subj', 'age', 'started'],
      yColumnNames: ['subj', 'age', 'started'],
    });
    await awaitCheck(() => document.querySelector('.d4-matrix-plot') !== null, 'matrix plot not found', 3000);

    if (matrixPlot.props.xColumnNames.length != 3)
      throw 'xColumnNames property was set incorrectly'; 
    if (matrixPlot.props.yColumnNames.length != 3)
      throw 'yColumnNames property was set incorrectly';     
  
    matrixPlot.setOptions({
      title: 'Test Matrix Plot',
      cellPlotType: 'scatter',
    });

    await awaitCheck(() => matrixPlot.props.title === 'Test Matrix Plot', 'title property has not been set', 1000);
    if (matrixPlot.props.cellPlotType != 'scatter')
      throw 'cellPlotType property has not been set to "scatter"'; 
  });

  // Does not work through Test Manager
  test('matrixPlot.serialization', async () => {
    await uploadProject('Test project with Matrix Plot', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Matrix Plot');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Matrix plot');
    const matrixPlot = findViewer('Matrix plot', v);
    
    if (matrixPlot!.props.xColumnNames.length != 3)
      throw 'xColumnNames property has not been deserialized incorrectly'; 
    if (matrixPlot!.props.yColumnNames.length != 3)
      throw 'yColumnNames property has not been deserialized incorrectly';
    if (matrixPlot!.props.title != 'Test Matrix Plot')
      throw 'title property has not been deserialized';
    if (matrixPlot!.props.cellPlotType != 'scatter')
      throw 'cellPlotType property has not been deserialized'; 
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Matrix Plot').first());
  }); 
});
