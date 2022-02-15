import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import {isColumnPresent, isViewerPresent, isExceptionElement, isErrorBallon} from './gui-utils';

category('Dialog: Cluster', () => {
    let v: DG.TableView;
    let demog = grok.data.demo.demog(1000);

    before(async () => {
        v = grok.shell.addTableView(demog);
    });

   test('dialogs.cluster', async () => {
        grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown'));
        grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown'));
        //isExceptionElement('Open Cluster Dialog (firs time)');
        await delay(1000);
        if (document.getElementsByClassName('d4-dialog').length == 0)
            throw 'Cluster dialog was not opened (first time)'

        let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);
        //isExceptionElement('Execute Cluster dialog (first time)');

        isColumnPresent(demog.columns, 'clusters');

        grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').root.dispatchEvent(new MouseEvent('mousedown'))
        //isExceptionElement('Open Cluster Dialog (second time');

        if (document.getElementsByClassName('d4-dialog').length == 0)
            throw 'Cluster dialog was not opened (second time)';

        let normalizeSelectorInput = document.getElementsByClassName('ui-input-choice ui-input-root')[1].childNodes[1] as HTMLSelectElement;
        normalizeSelectorInput.value = 'Z-scores';
        //isExceptionElement('Changed value of the Normalize field to Z-scores');

        let scatterPlotInput = document.getElementsByClassName('ui-input-bool ui-input-root')[0].childNodes[1] as HTMLElement;
        scatterPlotInput.click(); await delay(2000);
        //isExceptionElement('ScatterPlot field checked');

        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);
        //isExceptionElement('Execute Cluster dialog (second time)');

        isViewerPresent(Array.from(v.viewers), 'Scatter plot')
        isColumnPresent(demog.columns, 'clusters (2)');

   });

    /* after(async () => {
        v.close();
        grok.shell.closeAll();
    }); */

});