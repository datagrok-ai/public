import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from '../gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Viewers: Word Cloud', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('wordCloud.visual', async () => {
    let wordCloudIcon = document.getElementsByClassName('svg-word-cloud')[0] as HTMLElement;
    wordCloudIcon.click(); await delay(1000);

    isViewerPresent(Array.from(v.viewers), 'Word cloud');

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click(); await delay(1000);

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Word cloud') {
            throw 'Word Cloud viewer was not closed'
        }
    }
  }); 
  test('wordCloud.api', async () => {
    let wordCloud = v.addViewer(DG.VIEWER.WORD_CLOUD, {
        "wordColumnName": "study",
        "sizeColumnName": "study",
        "sizeColumnAggrType": "values"
      }); await delay(500);

    if (wordCloud.props.wordColumnName != 'study')
        throw 'Word column has not been set' 
    if (wordCloud.props.sizeColumnName != 'study')
        throw 'Size column has not been set'     
    if (wordCloud.props.sizeColumnAggrType != 'values')
        throw 'sizeColumnAggrType has not been set'           

    wordCloud.setOptions({
        title: 'Test Word Cloud',
        font: 'Time New Roman',
        maxWords: 30 
    }); await delay(500);

    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Word Cloud')
        throw 'title property has not been set' 

    if (wordCloud.props.font != 'Time New Roman')
        throw 'font property has not been set to Time New Roman' 
    if (wordCloud.props.maxWords != 30)
        throw 'maxWords property has not been set to 30'
  });  
  test('wordCloud.serialization', async () => {
    let project = DG.Project.create();
    project.name = 'Test project with Word Cloud'
    project.addChild(demog.getTableInfo());
    project.addChild(v.saveLayout());  
    await grok.dapi.layouts.save(v.saveLayout());
    await grok.dapi.tables.uploadDataFrame(demog);
    await grok.dapi.tables.save(demog.getTableInfo());
    await grok.dapi.projects.save(project);

    grok.shell.closeAll(); await delay(500);

    console.log('Project uploaded')

    await grok.dapi.projects.open('Test project with Word Cloud');

    console.log('after opening project')

    v = grok.shell.getTableView('demog 1000');

    isViewerPresent(Array.from(v.viewers), 'Word cloud');

    let wordCloud:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'Word cloud') {
            wordCloud = Array.from(v.viewers)[i];
            break;
        }
    }

    if (wordCloud!.props.wordColumnName != 'study')
        throw 'Word column has not been deserialized' 
    if (wordCloud!.props.sizeColumnName != 'study')
        throw 'Size column has not been deserialized'     
    if (wordCloud!.props.sizeColumnAggrType != 'values')
        throw 'sizeColumnAggrType has not been deserialized'
    if (wordCloud!.props.font != 'Time New Roman')
        throw 'font property has not been deserialized' 
    if (wordCloud!.props.maxWords != 30)
        throw 'maxWords property has not been deserialized'      
    let titleElem = document.querySelector("#elementContent > div.d4-layout-top > div > textarea") as HTMLSelectElement
    if (titleElem.value != 'Test Word Cloud')
        throw 'title property has not been deserialized' 
  }); 
  after(async () => {
    v.close();
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Word Cloud').first())
  }); 
});
