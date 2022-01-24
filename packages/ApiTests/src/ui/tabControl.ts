import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Tab control', () => {
    let v: DG.View;
    const tabs = ui.tabControl();
    tabs.addPane('First',()=>ui.div());
  
    before(async () => {
      v = grok.shell.newView('');
    });

    test('tabControl.root', async () => {
        checkHTMLElement('root', tabs.root, v, '.d4-tab-host')
    })

    test('tabControl.header', async () => {
        checkHTMLElement('header', tabs.header, v, '.d4-tab-header-stripe')
    })

    test('tabControl.pane', async () => {
        checkHTMLElement('header', tabs.header, v, '.d4-tab-header-stripe')
    })

    after(async () => {
        v.close();
        grok.shell.closeAll();
      });  

});