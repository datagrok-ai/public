import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import './gui-utils'
import {isDialogPresent, isViewerPresent, isExceptionElement, isErrorBallon, setDialogInputValue, returnDialog, isColumnPresent} from './gui-utils';
import {Column, ColumnList} from "datagrok-api/dg";

category('Dialog: Multivariate Analysis', () => {
    let v: DG.TableView;
    let demog = grok.data.demo.demog(1000);

    before(async () => {
        v = grok.shell.addTableView(demog);
    });

    test('dialogs.multivariateAnalysis', async () => {
        grok.shell.topMenu.find('Tools').find('Data Science').find('Multivariate Analysis (PLS)...').root.dispatchEvent(new MouseEvent('mousedown'));
        grok.shell.topMenu.find('Tools').find('Data Science').find('Multivariate Analysis (PLS)...').root.dispatchEvent(new MouseEvent('mousedown'));
        await delay(1000);
        isDialogPresent('Multivariate Analysis (PLS)');

        let okButton:HTMLElement | undefined;
        let button;
        for(let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++){
            button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
            if(button.innerText == 'OK') {
                okButton = button;
            }
        }
        if (okButton != undefined){
            if (okButton.classList[2] != 'disabled')
                throw 'OK Button does not have disabled class';
        }

        returnDialog('Multivariate Analysis (PLS)')!.input('Features').input.dispatchEvent(new MouseEvent('click')); await delay(1000);
        let okButtonInSelectColumns:HTMLElement | undefined;
        let allButtonInSelectColumns:HTMLElement | undefined;

        for(let i=0; i<document.getElementsByClassName('d4-link-label').length; i++) {
            button = document.getElementsByClassName('d4-link-label')[i] as HTMLElement;
            if (button.innerText == 'All') {
                allButtonInSelectColumns = button;
            }
        }

        for(let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++){
            button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
            if(button.innerText == 'OK') {
                okButtonInSelectColumns = button;
            }
        }
        if (allButtonInSelectColumns != undefined){
            allButtonInSelectColumns.click(); await delay(500);
        }
        if (okButtonInSelectColumns != undefined){
            okButtonInSelectColumns.click(); await delay(500);
        }
        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);
    });
    /* after(async () => {
        v.close();
        grok.shell.closeAll();
    }); */
});