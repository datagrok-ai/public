import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: PC Plot', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('pcPlot.visual', async () => {
    let pcPlotIcon = document.getElementsByClassName('svg-pc-plot')[0] as HTMLElement;
    pcPlotIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'PC Plot');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'PC Plot') {
            throw 'PC Plot viewer was not closed'
        }
    }
  }); 
  test('pcPlot.api', async () => {
    let pcPlot = v.addViewer(DG.VIEWER.PC_PLOT, {
        "columnNames": [
          "age",
          "sex",
          "started"
        ]
      }); await delay(500);

    if (pcPlot.props.columnNames.length != 3)
        throw 'columnNames property was set incorrectly' 
  
    pcPlot.setOptions({
        title: 'Test PC Plot',
        colorColumnName: "height",
        lineWidth: 3,
        showFilters: false
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test PC Plot')
        throw 'title property has not been set'

    if (pcPlot.props.colorColumnName != 'height')
        throw 'colorColumnName property has not been set to "height"' 
    if (pcPlot.props.lineWidth != 3)
        throw 'lineWidth property has not been set' 
    if (pcPlot.props.showFilters)
        throw 'showFilters property has not been set to "false"' 

    let rangeSelectors = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-pc-plot ui-box')[0].getElementsByTagName('svg');
    for (let i = 0; i < rangeSelectors.length; i++){
        if (rangeSelectors[i].style.visibility != 'hidden'){
            throw 'rangeSelectors not hidden after disabling "showFilters" property'
        }
    }         
  });  
  test('pcPlot.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with PC Plot'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with P C Plot');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'PC Plot');

    let pcPlot:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'PC Plot') {
            pcPlot = Array.from(v.viewers)[i];
            break;
        }
    }
    if (pcPlot!.props.columnNames.length != 3)
        throw 'columnNames property has not been deserialized' 
    if (pcPlot!.props.colorColumnName != 'height')
     throw 'colorColumnName property has not been deserialized' 
    if (pcPlot!.props.lineWidth != 3)
        throw 'lineWidth property has not been deserialized' 
    if (pcPlot!.props.showFilters)
        throw 'showFilters property has not been deserialized'
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with P C Plot').first())
  }); 
});
