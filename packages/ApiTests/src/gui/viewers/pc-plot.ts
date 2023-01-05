/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: PC Plot', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('pcPlot.visual', async () => {
    const pcPlotIcon = document.getElementsByClassName('svg-pc-plot')[0] as HTMLElement;
    pcPlotIcon.click();
    await awaitCheck(() => document.querySelector('.d4-pc-plot') !== null, 'pc plot not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'PC Plot');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click();
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null, 'pc plot properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'PC Plot viewer was not closed', 1000);
  });

  test('pcPlot.api', async () => {
    const pcPlot = v.addViewer(DG.VIEWER.PC_PLOT, {
      'columnNames': [
        'age',
        'sex',
        'started',
      ],
    });
    await awaitCheck(() => document.querySelector('.d4-pc-plot') !== null, 'pc plot not found', 3000);

    if (pcPlot.props.columnNames.length != 3)
      throw 'columnNames property was set incorrectly'; 
  
    pcPlot.setOptions({
      title: 'Test PC Plot',
      colorColumnName: 'height',
      lineWidth: 3,
      showFilters: false,
    });

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test PC Plot', 'title property has not been set', 2000);
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
    await uploadProject('Test project with P C Plot', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with P C Plot');
    v = grok.shell.getTableView('demog 1000');
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
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with P C Plot').first());
  }); 
});
