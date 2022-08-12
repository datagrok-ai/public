import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, uploadProject, findViewer} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Pie Chart', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('pieChart.visual', async () => {
    let pieChartIcon = document.getElementsByClassName('svg-pie-chart')[0] as HTMLElement;
    pieChartIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Pie chart');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'


    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Pie chart') {
            throw 'Pie chart viewer was not closed'
        }
    }
  }); 
  test('pieChart.api', async () => {
    let pieChart = v.pieChart({
        category: 'race',
    }); await delay(500);

    if (pieChart.props.categoryColumnName != 'race')
        throw 'Category column has not been set' 
    
    pieChart.setOptions({
        title: 'Test Pie chart',
        includeNulls: false,
        segmentLengthColumnName: 'age',
        onClick: 'Filter',
        startAngle: 30 
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Pie chart')
        throw 'title property has not been set'
    
    if (pieChart.props.includeNulls)
        throw 'includeNulls property has not been set to false' 
    if (pieChart.props.segmentLengthColumnName != 'age')
        throw 'segmentLengthColumnName property has not been set'
    if (pieChart.props.onClick != 'Filter')
        throw 'onClick property has not been set'
    if (pieChart.props.startAngle != 30)
        throw 'startAngle property (default) has not been set'
  });  
  test('pieChart.serialization', async () => {
    await uploadProject('Test project with Pie Chart', demog.getTableInfo(), v, demog);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Pie Chart');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Pie chart');

    let pieChart = findViewer('Pie chart', v);

    if (pieChart!.props.includeNulls)
        throw 'includeNulls property has not been deserialized' 
    if (pieChart!.props.segmentLengthColumnName != 'age')
        throw 'segmentLengthColumnName property has not been deserialized'
    if (pieChart!.props.onClick != 'Filter')
        throw 'onClick property has not been deserialized'
    if (pieChart!.props.startAngle != 30)
        throw 'startAngle property (default) has not been deserialized'
    if (pieChart!.props.categoryColumnName != 'race')
        throw 'Category column has not been deserialized' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Pie Chart').first())
  }); 
});
