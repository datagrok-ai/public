import {after, before, category, delay, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('UI - Inputs', () => {
    let v: DG.View;
    const inputs = [
        ui.stringInput('',''),
        ui.intInput('', 0),
        ui.floatInput('', 0.00),
        ui.boolInput('',true),
        ui.dateInput('', DG.DateTime.fromDate(new Date(2000, 1, 1))),
        ui.textInput('',''),
        ui.searchInput('','')
    ];

    before(async () => {
        v = grok.shell.newView('');
    });

    test('input.root', async () => {
        for (const input of inputs){
            checkHTMLElement(input, '.ui-input-root');
        }
    })

    test('input.input', async () => {
        for (const input of inputs){
            checkHTMLElement(input, '.ui-input-editor');
        }
    })

    test('input.captionLabel', async () => {
        for (const input of inputs){
            checkHTMLElement(input, '.ui-input-label');
        }
    })

    test('input.caption', async () => {
        for (const input of inputs){
            checkStringValue(input, '.ui-input- ');
        }
    })

    test('input.stringValue', async () => {
        for (const input of inputs){
            checkStringValue(input, '.ui-input-editor');
        }
    })



    after(async () => {
        v.close();
        grok.shell.closeAll();
    });

    function checkHTMLElement(input: DG.InputBase, selector: string): void {
        v.append(input.root);
        let e = v.root.querySelector(selector);
        if (e == undefined)
        throw `Element "${selector}" not found`;
        input.root.remove();
    }

    function checkStringValue(input: DG.InputBase, selector: string): void {
        v.append(input.root);
        let value = (<HTMLInputElement>v.root.querySelector(selector)).innerText;
        try {
            expect(input.caption, value)
        }
        catch (x){
            throw x;
        }
        finally{
            input.root.remove();
        }
    }

});

