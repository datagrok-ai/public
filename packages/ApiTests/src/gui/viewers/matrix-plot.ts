import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, uploadProject, findViewer} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Matrix Plot', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('matrixPlot.visual', async () => {
    let matrixPlotIcon = document.getElementsByClassName('svg-matrix-plot')[0] as HTMLElement;
    matrixPlotIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Matrix plot');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Matrix Plot') {
            throw 'Matrix Plot viewer was not closed'
        }
    }
  }); 
  test('matrixPlot.api', async () => {
    let matrixPlot = v.matrixPlot({
        xColumnNames: ['subj', 'age', 'started'],
        yColumnNames: ['subj', 'age', 'started']
        }); await delay(500);

    if (matrixPlot.props.xColumnNames.length != 3)
        throw 'xColumnNames property was set incorrectly' 
    if (matrixPlot.props.yColumnNames.length != 3)
        throw 'yColumnNames property was set incorrectly'     
  
    matrixPlot.setOptions({
        title: 'Test Matrix Plot',
        cellPlotType: 'scatter',
    }); await delay(500);

    if (matrixPlot.props.title != 'Test Matrix Plot')
        throw 'title property has not been set'
    if (matrixPlot.props.cellPlotType != 'scatter')
        throw 'cellPlotType property has not been set to "scatter"' 
  });  
  test('matrixPlot.serialization', async () => {
    await uploadProject('Test project with Matrix Plot', demog.getTableInfo(), v, demog);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Matrix Plot');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Matrix plot');

    let matrixPlot = findViewer('Matrix plot', v);
    
    if (matrixPlot!.props.xColumnNames.length != 3)
        throw 'xColumnNames property has not been deserialized incorrectly' 
    if (matrixPlot!.props.yColumnNames.length != 3)
        throw 'yColumnNames property has not been deserialized incorrectly'
        if (matrixPlot!.props.title != 'Test Matrix Plot')
        throw 'title property has not been deserialized'
    if (matrixPlot!.props.cellPlotType != 'scatter')
        throw 'cellPlotType property has not been deserialized' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Matrix Plot').first())
  }); 
});
