import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';

category('Dialog: Cluster', () => {
    let v: DG.TableView;
    let demog = grok.data.demo.demog(1000);

    before(async () => {
        v = grok.shell.addTableView(demog);
    });

   test('dialogs.cluster', async () => {
        //grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown'));
        grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Cluster Data');

        let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        isColumnPresent(demog.columns, 'clusters');

        grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Cluster Data');

        returnDialog('Cluster Data')!.input('Show scatter plot').input.click(); await delay(2000);

        isViewerPresent(Array.from(v.viewers), 'Scatter plot');
        isColumnPresent(demog.columns, 'clusters (2)');

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
             if (Array.from(v.viewers)[i].type == 'Scatter Plot'){
                throw 'Scatter Plot did not disappear after clicking on the "Cancel" button';
             }
        }

        if (demog.columns.byName('clusters (2)') != null)
            throw 'cluster (2) column did not disappear after clicking on the "Cancel" button';

        grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('Cluster Data');

        setDialogInputValue('Cluster Data', 'Normalize', 'Z-scores'); await delay(500);
        setDialogInputValue('Cluster Data', 'Clusters', 4); await delay(500);
        setDialogInputValue('Cluster Data', 'Metric', 'Manhattan'); await delay(500);
        returnDialog('Cluster Data')!.input('Show scatter plot').input.click(); await delay(2000);

        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        isViewerPresent(Array.from(v.viewers), 'Scatter plot')
        isColumnPresent(demog.columns, 'clusters (2)');

   });

     after(async () => {
        v.close();
        grok.shell.closeAll();
    });

});