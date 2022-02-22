import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import './gui-utils'
import {isDialogPresent, isViewerPresent, isExceptionElement, isErrorBallon, setDialogInputValue, returnDialog, isColumnPresent} from './gui-utils';
import {Column, ColumnList} from "datagrok-api/dg";

category('Dialog: Random Data', () => {
    let v: DG.TableView;
    let demog = grok.data.demo.demog(1000);

    before(async () => {
        v = grok.shell.addTableView(demog);
    });

    test('dialogs.randomData', async () => {
        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown'));
        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        setDialogInputValue('Random Data', 'Distribution', 'log-normal'); await delay(500);
        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        setDialogInputValue('Random Data', 'Distribution', 'binomial'); await delay(500);
        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        setDialogInputValue('Random Data', 'Distribution', 'poisson'); await delay(500);
        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        setDialogInputValue('Random Data', 'Distribution', 'uniform'); await delay(500);
        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(3000);

        isColumnPresent(demog.columns, 'normal');
        isColumnPresent(demog.columns, 'binomial');
        isColumnPresent(demog.columns, 'poisson');
        isColumnPresent(demog.columns, 'uniform');
        isColumnPresent(demog.columns, 'log-normal');

        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        returnDialog('Random Data')!.input('Show histogram').input.click(); await delay(3000);
        isViewerPresent(Array.from(v.viewers), 'Histogram');

        let cancelButton:HTMLElement;
        let button;
        for(let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++){
            button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
            if(button.innerText == 'CANCEL') {
                cancelButton = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
                cancelButton.click();
            }
        }

        await delay(1000);

        for(let i:number = 0; i < Array.from(v.viewers).length; i++){
            if (Array.from(v.viewers)[i].type == 'Histogram'){
                throw 'Histogram did not disappear after clicking on the "Cancel" button';
            }
        }

        grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Random Data');
        returnDialog('Random Data')!.input('Show histogram').input.click(); await delay(3000);
        isViewerPresent(Array.from(v.viewers), 'Histogram');
        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(3000);

        isViewerPresent(Array.from(v.viewers), 'Histogram')
    });

    after(async () => {
        v.close();
        grok.shell.closeAll();
    }); 
});