/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, settled, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Density plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
  });

  test('densityPlot.visual', async () => {
    v = grok.shell.addTableView(demog);
    const densityPlotIcon = document.getElementsByClassName('svg-density-plot')[0] as HTMLElement;
    densityPlotIcon.click(); 
    await waitUntil(() => document.querySelector('.d4-density-plot') !== null, 'density plot not found');
    isViewerPresent(Array.from(v.viewers), 'Density plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await waitUntil(() => document.querySelector('.grok-prop-panel') !== null,
      'density plot properties not found');
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); 
    await waitUntil(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found');
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await waitUntil(() => Array.from(v.viewers).length === 1, 'Density plot viewer was not closed');
  });

  test('densityPlot.api', async () => {
    v = grok.shell.addTableView(demog);
    const densityPlot = v.addViewer(DG.VIEWER.DENSITY_PLOT, {
      'xColumnName': 'weight',
      'yColumnName': 'height',
    }); 
    await waitUntil(() => document.querySelector('.d4-density-plot') !== null, 'density plot not found');

    if (densityPlot.props.xColumnName != 'weight')
      throw 'xColumnName property was set incorrectly'; 
    if (densityPlot.props.yColumnName != 'height')
      throw 'yColumnName property was set incorrectly'; 

    densityPlot.setOptions({
      title: 'Test Density plot',
      showXAxis: false,
      bins: 29,
    });

    await waitUntil(() => (Array.from(document.querySelectorAll(
      '#elementContent > div.d4-layout-top > div > textarea')) as HTMLTextAreaElement[])
      .some((e) => e.value === 'Test Density plot'), 'title property has not been set');
    if (densityPlot.props.showXAxis)
      throw 'showXAxis property has not been set to false'; 
    if (densityPlot.props.bins != 29)
      throw 'bins property has not been set'; 
  });

  // Does not work through Test Manager
  test('densityPlot.serialization', async () => {
    v = grok.shell.addTableView(demog);
    await settled(v.addViewer(DG.VIEWER.DENSITY_PLOT, {
      xColumnName: 'weight',
      yColumnName: 'height',
      title: 'Test Density plot',
      showXAxis: false,
      bins: 29,
    }));
    const project = await uploadProject('Test project with Density plot', demog.getTableInfo(), v, demog);
    try {
      grok.shell.closeAll();
      await project.open();
      v = grok.shell.tableView('demog 1000');
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
      if (densityPlot!.props.bins != 29)
        throw 'bins property has not been deserialized'; 
    } finally {
      await grok.dapi.projects.delete(project);
    }
  });
}, { owner: 'agolovko@datagrok.ai' });
