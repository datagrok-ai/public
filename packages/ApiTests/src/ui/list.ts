import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: List', () => {
    let v: DG.View;

    let list = ui.list([
        'element 1',
        'element 2'
    ]);


    before(async () => {
        v = grok.shell.newView('');
    });

    test('list.root', async () => {
        checkHTMLElement('list', list, v, '.d4-flex-col');
    });

    after(async () => {
        v.close();
        grok.shell.closeAll();
    });

});