import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Image', () => {
    let v: DG.View;

    let image = ui.image('https://datagrok.ai/help/visualize/viewers-interaction-main.gif',400,200);

    before(async () => {
        v = grok.shell.newView('');
    });

    test('image.root', async () => {
        checkHTMLElement('image', image, v, '.ui-image');
    });

    test('image.click', async () => {
        onClick(image);
    });

    test('image.size', async () => {
        if(image.style.width != '400px' || image.style.height != '200px')
            throw 'image size error'
    });

    after(async () => {
        v.close();
        grok.shell.closeAll();
    });

    function onClick(root: HTMLElement): void {
        v.append(root);

        let check = false;
        root.addEventListener("click", function () {
            check = true;
        });

        root.click();

        if (check == false)
            throw `"${root}": OnClick error`;

        root.remove();
    }

});