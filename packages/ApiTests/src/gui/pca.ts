import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from '../ui/utils';
import './gui-utils'
import {isDialogPresent, isViewerPresent, isExceptionElement, isErrorBallon, setDialogInputValue, returnDialog, isColumnPresent} from './gui-utils';
import {Column, ColumnList} from "datagrok-api/dg";

category('Dialog: PCA', () => {
    let v: DG.TableView;
    let demog = grok.data.demo.demog(1000);

    before(async () => {
        v = grok.shell.addTableView(demog);
    });

    test('dialogs.pca', async () => {
        grok.shell.topMenu.find('Tools').find('Data Science').find('Principal Component Analysis...').root.dispatchEvent(new MouseEvent('mousedown'));
        grok.shell.topMenu.find('Tools').find('Data Science').find('Principal Component Analysis...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('PCA');

        let okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);
        isErrorBallon('Errors calling PCA: features: Value not defined.');

        grok.shell.topMenu.find('Tools').find('Data Science').find('Principal Component Analysis...').root.dispatchEvent(new MouseEvent('mousedown')); await delay(1000);
        isDialogPresent('PCA');

        let featuresField:DG.Column[] = [];
        featuresField.push(demog.columns.byName('age'));
        featuresField.push(demog.columns.byName('height'));
        featuresField.push(demog.columns.byName('weight'));
        setDialogInputValue('PCA','Features',featuresField);

        //Open ColumnSelector for Features field:
        returnDialog('PCA')!.input('Features').input.dispatchEvent(new MouseEvent('click')); await delay(1000);
        let okButtonInSelectColumn = returnDialog('Select columns...')?.root.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButtonInSelectColumn.click(); await delay(1000)

        setDialogInputValue('PCA', 'Components', 3); await delay(500);
        returnDialog('PCA')!.input('Center').input.click(); await delay(500);

        okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButton!.click(); await delay(1000);

        isColumnPresent(demog.columns, 'PCA0');
        isColumnPresent(demog.columns, 'PCA1');
        isColumnPresent(demog.columns, 'PCA2');
    });

     after(async () => {
        v.close();
        grok.shell.closeAll();
    });
});