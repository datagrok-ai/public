import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Bar Chart', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('barChart.visual', async () => {
    let barChartIcon = document.getElementsByClassName('svg-bar-chart')[0] as HTMLElement;
    barChartIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Bar chart');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'


    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Bar chart') {
            throw 'Bar chart viewer was not closed'
        }
    }
  });
  test('barChart.api', async () => {
    let barChart = v.barChart({
        split: 'race',
        value: 'height',
        valueAggrType: 'max'
    }); await delay(1000);

    if (barChart.props.splitColumnName != 'race')
        throw 'Split column has not been set' 
    if (barChart.props.valueColumnName != 'height')
        throw 'Value column has not been set'
    if (barChart.props.valueAggrType != 'max')
        throw 'Aggregation has not been set' 
    
    barChart.setOptions({
        title: 'Test Bar chart',
        relativeValues: true,
        onClick: 'Filter',
        barBorderLineMouseOverWidth: 10 
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Bar chart')
        throw 'title property has not been set'
    
    if (!barChart.props.relativeValues)
        throw 'relativeValues property has not been set' 
    if (barChart.props.onClick != 'Filter')
        throw 'onClick property has not been set'
    if (barChart.props.barBorderLineMouseOverWidth != 10)
        throw 'barBorderLineMouseOverWidth property has not been set'
    if (barChart.props.barSortOrder != 'desc')
        throw 'barSortOrder property (default) has not been set'
    if (barChart.props.barSortType != 'by value')
        throw 'barSortType property (default) has not been set'   

  });
  test('barChart.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with Bar Chart'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Bar Chart');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Bar chart');

    let barChart:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Bar chart') {
            barChart = Array.from(v.viewers)[i];
            break;
        }
    }

    if (!barChart!.props.relativeValues)
        throw 'relativeValues property has not been deserialize' 
    if (barChart!.props.onClick != 'Filter')
        throw 'onClick property has not been deserialize'
    if (barChart!.props.barBorderLineMouseOverWidth != 10)
        throw 'barBorderLineMouseOverWidth property has not been deserialize'
    if (barChart!.props.barSortOrder != 'desc')
        throw 'barSortOrder property (default) has not been deserialize'
    if (barChart!.props.barSortType != 'by value')
        throw 'barSortType property (default) has not been deserialize'
    if (barChart!.props.splitColumnName != 'race')
        throw 'Split column has not been deserialize' 
    if (barChart!.props.valueColumnName != 'height')
        throw 'Value column has not been deserialize'
    if (barChart!.props.valueAggrType != 'max')
        throw 'Aggregation has not been deserialize' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
  }); 
});
