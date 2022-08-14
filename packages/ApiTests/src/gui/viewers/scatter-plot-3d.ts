import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, uploadProject, findViewer} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: 3D Scatter Plot', () => {  
  const demog = grok.data.demo.demog(1000);
  let v:DG.TableView;
  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('scatterPlot3D.visual', async () => {
    let ScatterPlot3DIcon = document.getElementsByClassName('svg-3d-scatter-plot')[0] as HTMLElement;
    ScatterPlot3DIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), '3d scatter plot');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == '3d scatter plot') {
            throw '3d scatter plot viewer was not closed'
        }
    }
  });
  test('scatterPlot3D.api', async () => {
    let scatterPlot3D = v.addViewer(DG.VIEWER.SCATTER_PLOT_3D, {
        "xColumnName": "height",
        "yColumnName": "weight",
        "zColumnName": "age",
        "sizeColumnName": "race",
        "colorColumnName": "sex"
      }); await delay(1000);

    if (scatterPlot3D.props.xColumnName != 'height')
        throw 'X column has not been set' 
    if (scatterPlot3D.props.yColumnName != 'weight')
        throw 'Y column has not been set'
    if (scatterPlot3D.props.zColumnName != 'age')
        throw 'Z column has not been set' 
    if (scatterPlot3D.props.sizeColumnName != 'race')
        throw 'Size column has not been set'
    if (scatterPlot3D.props.colorColumnName != 'sex')
        throw 'Color column has not been set' 
    
    scatterPlot3D.setOptions({
        title: 'Test 3D Scatter Plot',
        dynamicCameraMovement: true,
        markerType: "sphere",
        markerOpacity: 17
    }); await delay(500);

    if (scatterPlot3D.props.title != 'Test 3D Scatter Plot')
        throw 'title property has not been set'   
    if (!scatterPlot3D.props.dynamicCameraMovement)
        throw 'dynamicCameraMovement property has not been set' 
    if (scatterPlot3D.props.markerType != 'sphere')
        throw 'markerType property has not been set'
    if (scatterPlot3D.props.markerOpacity != 17)
        throw 'markerOpacity property has not been set'
  });
  test('scatterPlot3D.serialization', async () => {
    await uploadProject('Test project with 3 D Scatter Plot', demog.getTableInfo(), v, demog);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with 3 D Scatter Plot');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), '3d scatter plot');

    let scatterPlot3D = findViewer('3d scatter plot', v);

    if (!scatterPlot3D!.props.dynamicCameraMovement)
        throw 'dynamicCameraMovement property has not been deserialized' 
    if (scatterPlot3D!.props.markerType != 'sphere')
        throw 'markerType property has not been deserialized'
    if (scatterPlot3D!.props.markerOpacity != 17)
        throw 'markerOpacity property has not been deserialized'
    if (scatterPlot3D!.props.xColumnName != 'height')
        throw 'X column has not been deserialized' 
    if (scatterPlot3D!.props.yColumnName != 'weight')
        throw 'Y column has not been deserialized'
    if (scatterPlot3D!.props.zColumnName != 'age')
        throw 'Z column has not been deserialized' 
    if (scatterPlot3D!.props.sizeColumnName != 'race')
        throw 'Size column has not been deserialized'
    if (scatterPlot3D!.props.colorColumnName != 'sex')
        throw 'Color column has not been deserialized' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with 3 D Scatter Plot').first())
  }); 
});
