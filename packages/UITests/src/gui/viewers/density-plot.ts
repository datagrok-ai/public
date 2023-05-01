/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Density plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('densityPlot.visual', async () => {
    const densityPlotIcon = document.getElementsByClassName('svg-density-plot')[0] as HTMLElement;
    densityPlotIcon.click(); 
    await awaitCheck(() => document.querySelector('.d4-density-plot') !== null, 'density plot not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Density plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null,
      'density plot properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); 
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Density plot viewer was not closed', 1000);
  });

  test('densityPlot.api', async () => {
    const densityPlot = v.addViewer(DG.VIEWER.DENSITY_PLOT, {
      'xColumnName': 'weight',
      'yColumnName': 'height',
    }); 
    await awaitCheck(() => document.querySelector('.d4-density-plot') !== null, 'density plot not found', 3000);

    if (densityPlot.props.xColumnName != 'weight')
      throw 'xColumnName property was set incorrectly'; 
    if (densityPlot.props.yColumnName != 'height')
      throw 'yColumnName property was set incorrectly'; 

    densityPlot.setOptions({
      title: 'Test Density plot',
      showXAxis: false,
      xBins: 29,
    });

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Density plot', 'title property has not been set', 2000);
    if (densityPlot.props.showXAxis)
      throw 'showXAxis property has not been set to false'; 
    if (densityPlot.props.xBins != 29)
      throw 'xBins property has not been set'; 
  });

  // Does not work through Test Manager
  test('densityPlot.serialization', async () => {
    await uploadProject('Test project with Density plot', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Density plot');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Density plot');
    const densityPlot = findViewer('Density plot', v);

    if (densityPlot!.props.xColumnName != 'weight')
      throw 'xColumnName property has not been deserialized'; 
    if (densityPlot!.props.yColumnName != 'height')
      throw 'yColumnName property has not been deserialized'; 
    if (densityPlot!.props.title != 'Test Density plot')
      throw 'title property has not been deserialized';
    if (densityPlot!.props.showXAxis)
      throw 'showXAxis property has not been deserialized'; 
    if (densityPlot!.props.xBins != 29)
      throw 'xBins property has not been deserializedt'; 
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Density plot').first());
  }); 
});
