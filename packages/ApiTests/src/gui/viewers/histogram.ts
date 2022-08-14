import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, uploadProject, findViewer} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Histogram', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('histogram.visual', async () => {
    let histogramChartIcon = document.getElementsByClassName('svg-histogram')[0] as HTMLElement;
    histogramChartIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Histogram');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Histogram') {
            throw 'Histogram viewer was not closed'
        }
    }
  }); 
  test('histogram.api', async () => {
    let histogram = v.histogram({
        value: 'weight',
        binWidthRatio: 2
        }); await delay(500);

    if (histogram.props.valueColumnName != 'weight')
        throw 'Value column has not been set' 
    if (histogram.props.binWidthRatio != 2)
        throw 'binWidthRatio has not been set to 2'     

    let filterCheckbox = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0].getElementsByClassName('d4-filter-out-missing-values')[0] as HTMLInputElement;
    if (!filterCheckbox)
        throw 'filterCheckbox element not found'      

    histogram.setOptions({
        title: 'Test Histogram',
        colorColumnName: 'age',
        showRangeInputs: true,
        bins: 30 
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Histogram')
        throw 'title property has not been set' 

    if (histogram.props.bins != 30)
        throw 'bins property has not been set to 30' 
    if (histogram.props.colorColumnName != 'age')
        throw 'colorColumnName property has not been set'

    let rangeInputMin = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0].getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-min')[0] as HTMLInputElement;
    let rangeInputMax = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0].getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-max')[0] as HTMLInputElement;
    if (!rangeInputMin)
       throw 'min range input from showRangeInputs property not found'
    if (!rangeInputMax)
       throw 'max range input from showRangeInputs property not found'   

  });  
  test('histogram.serialization', async () => {
    await uploadProject('Test project with Histogram', demog.getTableInfo(), v, demog);

    grok.shell.closeAll(); await delay(500);

    await grok.dapi.projects.open('Test project with Histogram');

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Histogram');

    let histogram = findViewer('Histogramt', v);

    if (histogram!.props.valueColumnName != 'weight')
        throw 'Value column has not been deserialized' 
    if (histogram!.props.binWidthRatio != 2)
        throw 'binWidthRatio has not been deserialized'   
    if (histogram!.props.bins != 30)
        throw 'bins property has not been deserialized' 
    if (histogram!.props.colorColumnName != 'age')
        throw 'colorColumnName property has not been deserialized'
    if (!histogram!.props.showRangeInputs)
        throw 'showRangeInputs property has not been deserialized'    

    let filterCheckbox = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0].getElementsByClassName('d4-filter-out-missing-values')[0] as HTMLInputElement;
    if (!filterCheckbox)
        throw 'filterCheckbox element not found'  

    let rangeInputMin = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0].getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-min')[0] as HTMLInputElement;
    let rangeInputMax = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-histogram ui-box')[0].getElementsByClassName('ui-input-editor d4-filter-input d4-filter-input-max')[0] as HTMLInputElement;
    if (!rangeInputMin)
        throw 'min range input from showRangeInputs property not found'
    if (!rangeInputMax)
        throw 'max range input from showRangeInputs property not found' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Histogram').first())
  }); 
});
