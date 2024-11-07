import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from './utils';

category('UI: Tab control', () => {
  let v: DG.View;
  let tabs: DG.TabControl;

  before(async () => {
    const tabsItems = {
      '1' : ui.div([], 'new-pane'),
      '2' : ui.div([], 'new-pane')
    };
    v = grok.shell.newView('');
    tabs = ui.tabControl(tabsItems);
    v.append(tabs.root);
  });

  test('tabControl.root', async () => {
    checkHTMLElement('root', tabs.root, v, '.d4-tab-host');
  });

  test('tabControl.header', async () => {
    checkHTMLElement('header', tabs.header, v, '.d4-tab-header-stripe');
  });

  test('tabControl.addPane', async () => {
    tabs.addPane('New pane', ()=> ui.div([], 'new-pane'));
    if (tabs.panes.length == 0)
      throw new Error('addPane error');
  });

  test('tabControl.currentPane', async () => {
    const pane = tabs.panes[0];
    const currentPane = tabs.currentPane;
    expect(pane, currentPane);
  });

  test('pane.content', async () => {
    checkHTMLElement('header', tabs.panes[0].content, v, '.new-pane');
  });

  test('pane.header', async () => {
    checkHTMLElement('header', tabs.panes[0].header, v, '.d4-tab-header.selected');
  });

  test('tabControl.getPane', async () => {
    tabs.addPane('New pane', ()=> ui.div([], 'new-pane'));
    const pane = tabs.panes[tabs.panes.length-1];
    const getPane = tabs.getPane('New pane');
    expect(pane, getPane);
  });

  test('tabControl.onTabAdded', async () => {
    let check = true;
    tabs.onTabAdded.subscribe((_)=> {
      check=false;
    });
    tabs.addPane('New pane 2', ()=> ui.div([]));
    if (check)
      throw new Error('onTabAdded event error');
  });

  test('tabControl.onTabChanged', async () => {
    let check = true;
    tabs.onTabChanged.subscribe((_)=> {
      check=false;
    });
    tabs.currentPane = tabs.panes[1];
    if (check)
      throw new Error('onTabChanged event error');
  });

  test('tabControl.onBeforeTabChanged', async () => {
    tabs.currentPane = tabs.panes[0];
    let check = true;
    tabs.onBeforeTabChanged.subscribe((_)=> {check=false;});
    tabs.currentPane = tabs.panes[1];
    expect(check, false, 'onBeforeTabChanged event error')
  });

  after(async () => {
    grok.shell.closeAll();
  });
}, {clear: false});
