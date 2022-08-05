import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import './gui-utils';
import {isDialogPresent, isViewerPresent, isExceptionElement, isErrorBallon, setDialogInputValue, returnDialog, isColumnPresent} from './gui-utils';
import {Column, ColumnList} from 'datagrok-api/dg';

category('Dialog: Missing Values Imputation', () => {
  let v:DG.TableView;
  const demog = grok.data.loadTable('https://dev.datagrok.ai/demo/demog.csv');

  before(async () => {
    grok.shell.addTableView(await demog);
  });

  test('dialogs.missingValuesImputation', async () => {
    grok.shell.topMenu.find('Tools').find('Data Science').find('Missing Values Imputation...').click();
    await delay(1000);
    isDialogPresent('Missing Values Imputation');

    const imputeField:DG.Column[] = [];
    imputeField.push((await demog).columns.byName('AGE'));
    imputeField.push((await demog).columns.byName('HEIGHT'));
    imputeField.push((await demog).columns.byName('WEIGHT'));
    setDialogInputValue('Missing Values Imputation', 'Impute', imputeField);
        returnDialog('Missing Values Imputation')!.input('Impute').input.dispatchEvent(new MouseEvent('click')); await delay(500);
        let okButtonInSelectColumn = returnDialog('Select columns...')?.root.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButtonInSelectColumn.click(); await delay(500);

        const dataField:DG.Column[] = [];
        dataField.push((await demog).columns.byName('SEX'));
        dataField.push((await demog).columns.byName('RACE'));
        dataField.push((await demog).columns.byName('DIS_POP'));
        dataField.push((await demog).columns.byName('DEMOG'));
        setDialogInputValue('Missing Values Imputation', 'Data', dataField);
        returnDialog('Missing Values Imputation')!.input('Data').input.dispatchEvent(new MouseEvent('click')); await delay(500);
        okButtonInSelectColumn = returnDialog('Select columns...')?.root.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
        okButtonInSelectColumn.click(); await delay(500);

        let okButton:HTMLElement | undefined;
        let button;
        for (let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++) {
          button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
          if (button.innerText == 'OK')
            okButton = button;
        }
        if (okButton) {
          okButton.click();
          await delay(8000);
        }

        if (((await demog).col('AGE')?.stats.missingValueCount) != 0)
          throw 'column AGE contains empty rows!';

        if (((await demog).col('WEIGHT')?.stats.missingValueCount) != 0)
          throw 'column WEIGHT contains empty rows!';
        if (((await demog).col('HEIGHT')?.stats.missingValueCount) != 0)
          throw 'column HEIGHT contains empty rows!';
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
