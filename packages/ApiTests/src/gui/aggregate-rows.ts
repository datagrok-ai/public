import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';

category('Dialog: Aggregate Rows', () => {
    let v: DG.TableView;
    let demog = grok.data.demo.demog(1000);

    before(async () => {
        v = grok.shell.addTableView(demog);
    });

   test('dialogs.cluster', async () => {
        grok.shell.topMenu.find('Tools').find('Data').find('Aggregate Rows...').root.click(); await delay(1000);
        //isDialogPresent('Cluster Data');

        //let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        //okButton!.click(); await delay(1000);


   });

     after(async () => {
        v.close();
        grok.shell.closeAll();
    });

});