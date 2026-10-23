/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, settled, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: PC Plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
  });

  test('pcPlot.visual', async () => {
    v = grok.shell.addTableView(demog);
    const pcPlotIcon = document.getElementsByClassName('svg-pc-plot')[0] as HTMLElement;
    pcPlotIcon.click();
    await waitUntil(() => document.querySelector('.d4-pc-plot') !== null, 'pc plot not found');
    isViewerPresent(Array.from(v.viewers), 'PC Plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await waitUntil(() => document.querySelector('.grok-prop-panel') !== null, 'pc plot properties not found');
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await waitUntil(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found');
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await waitUntil(() => Array.from(v.viewers).length === 1, 'PC Plot viewer was not closed');
  });

  test('pcPlot.api', async () => {
    v = grok.shell.addTableView(demog);
    const pcPlot = v.addViewer(DG.VIEWER.PC_PLOT, {
      'columnNames': [
        'age',
        'sex',
        'started',
      ],
    });
    await waitUntil(() => document.querySelector('.d4-pc-plot') !== null, 'pc plot not found');

    if (pcPlot.props.columnNames.length != 3)
      throw 'columnNames property was set incorrectly'; 
  
    pcPlot.setOptions({
      title: 'Test PC Plot',
      colorColumnName: 'height',
      lineWidth: 3,
      showFilters: false,
    });

    await waitUntil(() => (Array.from(document.querySelectorAll(
      '#elementContent > div.d4-layout-top > div > textarea')) as HTMLTextAreaElement[])
      .some((e) => e.value === 'Test PC Plot'), 'title property has not been set');
    if (pcPlot.props.colorColumnName != 'height')
      throw 'colorColumnName property has not been set to "height"'; 
    if (pcPlot.props.lineWidth != 3)
      throw 'lineWidth property has not been set'; 
    if (pcPlot.props.showFilters)
      throw 'showFilters property has not been set to "false"'; 
    const rangeSelectors = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-pc-plot ui-box')[0]
      .getElementsByTagName('svg');
    for (let i = 0; i < rangeSelectors.length; i++) {
      if (rangeSelectors[i].style.visibility != 'hidden') 
        throw 'rangeSelectors not hidden after disabling "showFilters" property';
    }         
  });

  // Does not work through Test Manager
  test('pcPlot.serialization', async () => {
    v = grok.shell.addTableView(demog);
    await settled(v.addViewer(DG.VIEWER.PC_PLOT, {
      columnNames: ['age', 'sex', 'started'],
      title: 'Test PC Plot',
      colorColumnName: 'height',
      lineWidth: 3,
      showFilters: false,
    }));
    const project = await uploadProject('Test project with P C Plot', demog.getTableInfo(), v, demog);
    try {
      grok.shell.closeAll();
      await project.open();
      v = grok.shell.tableView('demog 1000');
      isViewerPresent(Array.from(v.viewers), 'PC Plot');
      const pcPlot = findViewer('PC Plot', v);

      if (pcPlot!.props.columnNames.length != 3)
        throw 'columnNames property has not been deserialized'; 
      if (pcPlot!.props.colorColumnName != 'height')
        throw 'colorColumnName property has not been deserialized'; 
      if (pcPlot!.props.lineWidth != 3)
        throw 'lineWidth property has not been deserialized'; 
      if (pcPlot!.props.showFilters)
        throw 'showFilters property has not been deserialized';
    } finally {
      await grok.dapi.projects.delete(project);
    }
  });
}, { owner: 'dkovalyov@datagrok.ai' });
