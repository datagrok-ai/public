import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Accordion', () => {
    let v: DG.View;
    const acc = ui.accordion();
    acc.addPane('Pane',()=>ui.div());
  
    before(async () => {
      v = grok.shell.newView('');
    });

    test('accordion.root', async () => {
        checkHTMLElement('accordion', acc.root, v, '.d4-accordion')
    })

    test('accordion-pane.root', async () => {
        checkHTMLElement('accordion pane', acc.root, v, '.d4-accordion-pane')
    })

    test('accordion-pane.name', async () => {
        if (acc.panes[0].name != 'Pane')
            throw 'pane.name error'
    })

    test('accordion-pane.expanded', async () => {
        if (acc.panes[0].expanded != false)
            throw 'pane.name error'
    })


    after(async () => {
        v.close();
        grok.shell.closeAll();
      });  

});