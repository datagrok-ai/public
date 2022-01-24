import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Accordion', () => {
    let v: DG.View;
    const acc = ui.accordion();
  
    before(async () => {
      v = grok.shell.newView('');
    });

    test('accordion.root', async () => {
        checkHTMLElement('accordion', acc.root, v, '.d4-accordion')
    })

    test('accordion.addPane', async () => {
        acc.addPane('New pane', ()=>{}, false);
        checkHTMLElement('accordion pane', acc.root, v, '.d4-accordion-pane')
    })

    test('accordion.getPane', async () => {
        if (acc.getPane('New pane') == undefined)
            throw 'getPane error'
    })

    test('pane.name', async () => {
        if (acc.panes[0].name != 'New pane')
            throw 'pane.name error'
    })

    test('pane.root', async () => {
        checkHTMLElement('Accordion Pane root', acc.panes[0].root, v, '.d4-accordion-pane')
    })

    test('pane.expanded', async () => {
        if (acc.panes[0].expanded != false)
            throw 'pane.name error'
    })

    test('accordion.removePane', async () => {
        acc.removePane(acc.panes[0])
        if (acc.getPane('New pane') != undefined)
            throw 'getPane error'
    })



    after(async () => {
        v.close();
        grok.shell.closeAll();
      });  

});