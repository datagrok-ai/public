import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Tab control', () => {
    let v: DG.View;
    let tabs = ui.tabControl();
    
    before(async () => {
      v = grok.shell.newView('');
    });

    test('tabControl.root', async () => {
        checkHTMLElement('root', tabs.root, v, '.d4-tab-host')
    })

    test('tabControl.header', async () => {
        checkHTMLElement('header', tabs.header, v, '.d4-tab-header-stripe')
    })

    test('tabControl.addPane', async () => {
        tabs.addPane('New pane', ()=> ui.div([],'new-pane'));
        if (tabs.panes.length == 0)
            throw 'addPane error'
    })

    test('tabControl.currentPane', async () => {
        let pane = tabs.panes[0];
        let currentPane = tabs.currentPane;
        expect(pane, currentPane)
    })

    test('pane.content', async () => {
        checkHTMLElement('header', tabs.panes[0].content, v, '.new-pane')
    })

    test('pane.header', async () => {
        checkHTMLElement('header', tabs.panes[0].header, v, '.d4-tab-header.selected')
    })

    test('tabControl.getPane', async () => {
        let pane = tabs.panes[0];
        let getPane = tabs.getPane('New pane');
        expect(pane, getPane)
    })

    test('tabControl.onTabAdded', async () => {
        let check = true;
        tabs.onTabAdded.subscribe((_)=> {
            check=false;
        });
        tabs.addPane('New pane 2', ()=> ui.div([]));
        if (check)
            throw 'onTabAdded event error'
    })

    test('tabControl.onTabChanged', async () => {
        let check = true;
        tabs.onTabChanged.subscribe((_)=> {
            check=false;
        });
        tabs.currentPane = tabs.panes[1]
        if (check)
            throw 'onTabChanged event error'
    })

    test('tabControl.onBeforeTabChanged', async () => {
        let check = true;
        tabs.onBeforeTabChanged.subscribe((_)=> {check=false;});
        tabs.currentPane = tabs.panes[0]
        if (check)
            throw 'onBeforeTabChanged event error'
    });

    after(async () => {
        v.close();
        tabs = ui.tabControl();
        grok.shell.closeAll();
      });  

});