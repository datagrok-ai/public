import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Box Plot', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('boxPlot.visual', async () => {
    let boxPlotIcon = document.getElementsByClassName('svg-box-plot')[0] as HTMLElement;
    boxPlotIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Box plot');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Box plot') {
            throw 'Box plot viewer was not closed'
        }
    }
  }); 
  test('boxPlot.api', async () => {
    let boxPlot = v.addViewer(DG.VIEWER.BOX_PLOT, {
        "valueColumnName": "age",
        "categoryColumnName": "race",
        "binColorColumnName": "height",
        "markerColorColumnName": "sex"
      }); await delay(500);

    if (boxPlot.props.valueColumnName != 'age')
        throw 'Value column has not been set' 
    if (boxPlot.props.categoryColumnName != 'race')
        throw 'Category column has not been set'     
    if (boxPlot.props.binColorColumnName != 'height')
        throw 'Bin color column has not been set'   
    if (boxPlot.props.markerColorColumnName != 'sex')
        throw 'Marker color column has not been set'    

    boxPlot.setOptions({
        title: 'Test Box Plot',
        markerSize: 10,
        showStatistics: true,
        axisType: 'logarithmic'
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Box Plot')
        throw 'title property has not been set' 

    if (boxPlot.props.markerSize != 10)
        throw 'Marker size property has not been set to 10' 
    if (!boxPlot.props.showStatistics)
        throw 'Show Statistics property has not been set to TRUE' 
    if (boxPlot.props.axisType != 'logarithmic')
        throw 'Axis Type property has not been set to logarithmic'
  });  
  test('boxPlot.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with Box Plot'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Box Plot');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Box plot');

    let boxPlot:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Box plot') {
            boxPlot = Array.from(v.viewers)[i];
            break;
        }
    }

    if (boxPlot!.props.valueColumnName != 'age')
        throw 'Value column has not been deserialized' 
    if (boxPlot!.props.categoryColumnName != 'race')
        throw 'Category column has not been deserialized'     
    if (boxPlot!.props.binColorColumnName != 'height')
        throw 'Bin color column has not been deserialized'   
    if (boxPlot!.props.markerColorColumnName != 'sex')
        throw 'Marker color column has not been deserialized'   
    if (boxPlot!.props.markerSize != 10)
        throw 'Marker size property has not been deserialized' 
    if (!boxPlot!.props.showStatistics)
        throw 'Show Statistics property hhas not been deserialized' 
    if (boxPlot!.props.axisType != 'logarithmic')
        throw 'Axis Type property has not been deserialized'
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Box Plot').first())
  }); 
});
