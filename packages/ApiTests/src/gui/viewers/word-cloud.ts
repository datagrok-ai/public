/* eslint-disable no-throw-literal */
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {isViewerPresent, uploadProject, findViewer} from '../gui-utils';


category('Viewers: Word Cloud', () => {
  let v: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
  });

  test('wordCloud.visual', async () => {
    const wordCloudIcon = document.getElementsByClassName('svg-word-cloud')[0] as HTMLElement;
    wordCloudIcon.click();
    await awaitCheck(() => document.querySelector('.d4-word-cloud') !== null, 'word cloud not found', 3000);
    isViewerPresent(Array.from(v.viewers), 'Word cloud');
    const protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); 
    await awaitCheck(() => document.querySelector('.grok-prop-panel') !== null,
      'word cloud properties not found', 1000);
    const hamburgerBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-menu')[0] as HTMLElement;
    hamburgerBtn.click();
    await awaitCheck(() => document.querySelector('.d4-menu-item-container.d4-vert-menu.d4-menu-popup') !== null,
      'hamburger menu not found', 1000);
    const closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0]
      .getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click();
    await awaitCheck(() => Array.from(v.viewers).length === 1, 'Word cloud viewer was not closed', 1000);
  });

  test('wordCloud.api', async () => {
    const wordCloud = v.addViewer(DG.VIEWER.WORD_CLOUD, {
      'wordColumnName': 'study',
      'sizeColumnName': 'study',
      'sizeColumnAggrType': 'values',
    });
    await awaitCheck(() => document.querySelector('.d4-word-cloud') !== null, 'word cloud not found', 3000);

    if (wordCloud.props.wordColumnName != 'study')
      throw 'Word column has not been set'; 
    if (wordCloud.props.sizeColumnName != 'study')
      throw 'Size column has not been set';     
    if (wordCloud.props.sizeColumnAggrType != 'values')
      throw 'sizeColumnAggrType has not been set';           

    wordCloud.setOptions({
      title: 'Test Word Cloud',
      font: 'Time New Roman',
      maxWords: 30, 
    });

    await awaitCheck(() => (document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement).
      value === 'Test Word Cloud', 'title property has not been set', 2000);
    if (wordCloud.props.font != 'Time New Roman')
      throw 'font property has not been set to Time New Roman'; 
    if (wordCloud.props.maxWords != 30)
      throw 'maxWords property has not been set to 30';
  });  

  // Does not work through Test Manager
  test('wordCloud.serialization', async () => {
    await uploadProject('Test project with Word Cloud', demog.getTableInfo(), v, demog);
    grok.shell.closeAll();
    await grok.dapi.projects.open('Test project with Word Cloud');
    v = grok.shell.getTableView('demog 1000');
    isViewerPresent(Array.from(v.viewers), 'Word cloud');
    const wordCloud = findViewer('Word cloud', v);

    if (wordCloud!.props.wordColumnName != 'study')
      throw 'Word column has not been deserialized'; 
    if (wordCloud!.props.sizeColumnName != 'study')
      throw 'Size column has not been deserialized';     
    if (wordCloud!.props.sizeColumnAggrType != 'values')
      throw 'sizeColumnAggrType has not been deserialized';
    if (wordCloud!.props.font != 'Time New Roman')
      throw 'font property has not been deserialized'; 
    if (wordCloud!.props.maxWords != 30)
      throw 'maxWords property has not been deserialized';      
    const titleElem = document.
      querySelector('#elementContent > div.d4-layout-top > div > textarea') as HTMLSelectElement;
    if (titleElem.value != 'Test Word Cloud')
      throw 'title property has not been deserialized'; 
  }); 

  after(async () => {
    grok.shell.closeAll();
    await grok.dapi.projects.delete(await grok.dapi.projects.filter('Test project with Word Cloud').first());
  }); 
});
