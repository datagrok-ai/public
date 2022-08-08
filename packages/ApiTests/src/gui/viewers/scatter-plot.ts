import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Scatter Plot', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('scatterPlot.visual', async () => {
    let scatterPlotIcon = document.getElementsByClassName('svg-scatter-plot')[0] as HTMLElement;
    scatterPlotIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Scatter plot');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);
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
    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Scatter plot') {
            throw 'viewer was not closed'
        }
    }
  });
  test('scatterPlot.api', async () => {
    let scatterPlot = v.scatterPlot({
        x: 'weight',
        y: 'height',
        size: 'age',
        color: 'sex'
    }); await delay(1000);

    if (scatterPlot.props.xColumnName != 'weight')
        throw 'X column has not been set' 
    if (scatterPlot.props.yColumnName != 'height')
        throw 'Y column has not been set'
    if (scatterPlot.props.sizeColumnName != 'age')
        throw 'Size column has not been set' 
    if (scatterPlot.props.colorColumnName != 'sex')
        throw 'Name column has not been set'

    scatterPlot.setOptions({
        title: 'Test Scatter Plot',
        showRegressionLine: true,
        markerType: 'dot',
        showVerticalGridLines: true,
        markerBorderWidth: 10 
    });

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Scatter Plot')
        throw 'title property has not been set'
    
    if (!scatterPlot.props.showRegressionLine)
        throw 'showRegressionLine property has not been set' 
    if (scatterPlot.props.markerType != 'dot')
        throw 'markerType property has not been set'
    if (!scatterPlot.props.showVerticalGridLines)
        throw 'showVerticalGridLines property has not been set'
    if (scatterPlot.props.markerBorderWidth != 10)
        throw 'markerBorderWidth property has not been set'
  });
  test('scatterPlot.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with Scatter plot'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Scatter plot');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Scatter plot');

    let scatterPlot:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Scatter plot') {
          scatterPlot = Array.from(v.viewers)[i];
          break;
        }
    }

    if (scatterPlot!.props.xColumnName != 'weight')
        throw 'X column has not been deserialized' 
    if (scatterPlot!.props.yColumnName != 'height')
        throw 'Y column has not been deserialized'
    if (scatterPlot!.props.sizeColumnName != 'age')
        throw 'Size column has not been deserialized' 
    if (scatterPlot!.props.colorColumnName != 'sex')
        throw 'Name column has not been deserialized'
    if (!scatterPlot!.props.showRegressionLine)
        throw 'showRegressionLine property has not been deserialized' 
    if (scatterPlot!.props.markerType != 'dot')
        throw 'markerType property has not been deserialized'
    if (!scatterPlot!.props.showVerticalGridLines)
        throw 'showVerticalGridLines property has not been deserialized'
    if (scatterPlot!.props.markerBorderWidth != 10)
        throw 'markerBorderWidth property has not been deserialized'         
  });
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Scatter plot').first())
  }); 
});
