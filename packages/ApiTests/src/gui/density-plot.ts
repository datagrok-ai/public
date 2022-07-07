import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Density plot', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });
  test('densityPlot.visual', async () => {
    let densityPlotIcon = document.getElementsByClassName('svg-density-plot')[0] as HTMLElement;
    densityPlotIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Density plot');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Density plot') {
            throw 'Density plot viewer was not closed'
        }
    }
  }); 
  test('densityPlot.api', async () => {
    let densityPlot = v.addViewer(DG.VIEWER.DENSITY_PLOT, {
        "xColumnName": "weight",
        "yColumnName": "height"
      }); await delay(500);

    if (densityPlot.props.xColumnName != 'weight')
        throw 'xColumnName property was set incorrectly' 
    if (densityPlot.props.yColumnName != 'height')
        throw 'yColumnName property was set incorrectly' 

    densityPlot.setOptions({
        title: 'Test Density plot',
        showXAxis: false,
        xBins: 29
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Density plot')
        throw 'title property has not been set'

    if (densityPlot.props.showXAxis)
        throw 'showXAxis property has not been set to false' 
    if (densityPlot.props.xBins != 29)
        throw 'xBins property has not been set' 
  });  
  test('densityPlot.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with Density plot'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Density plot');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Density plot');

    let densityPlot:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Density plot') {
            densityPlot = Array.from(v.viewers)[i];
            break;
        }
    }
    if (densityPlot!.props.xColumnName != 'weight')
        throw 'xColumnName property has not been deserialized' 
    if (densityPlot!.props.yColumnName != 'height')
        throw 'yColumnName property has not been deserialized' 
    if (densityPlot!.props.title != 'Test Density plot')
        throw 'title property has not been deserialized'
    if (densityPlot!.props.showXAxis)
        throw 'showXAxis property has not been deserialized' 
    if (densityPlot!.props.xBins != 29)
        throw 'xBins property has not been deserializedt' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Density plot').first())
  }); 
});
