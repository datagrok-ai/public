import {before, after, awaitCheck, category, test, delay} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog,
  setDialogInputValue, checkDialog} from './gui-utils';


category('GUI: Dialogs', () => {
  let v: DG.TableView;
  let df: DG.DataFrame;
  let demog: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(100);
    grok.shell.windows.simpleMode = true;
  });

  test('missingValuesImputation', async () => {
    demog = df.clone();
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.topMenu.find('Tools').find('Data Science').find('Missing Values Imputation...').click();
    await awaitCheck(() => checkDialog('Missing Values Imputation'), 'Dialog is not open', 1000);
    isDialogPresent('Missing Values Imputation');
    const imputeField: DG.Column[] = [];
    imputeField.push(demog.col('AGE') as DG.Column);
    imputeField.push(demog.col('HEIGHT') as DG.Column);
    imputeField.push(demog.col('WEIGHT') as DG.Column);
    setDialogInputValue('Missing Values Imputation', 'Impute', imputeField);
    returnDialog('Missing Values Imputation')!.input('Impute').input.dispatchEvent(new MouseEvent('click'));
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length === 2, 'cannot find Impute columns dialog', 1000);
    let okButtonInSelectColumn = returnDialog('Select columns...')?.root
      .getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButtonInSelectColumn.click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length === 1, 'error while selectiong impute columns', 1000);
    const dataField: DG.Column[] = [];
    dataField.push(demog.col('SEX') as DG.Column);
    dataField.push(demog.col('RACE') as DG.Column);
    dataField.push(demog.col('DIS_POP') as DG.Column);
    dataField.push(demog.col('DEMOG') as DG.Column);
    setDialogInputValue('Missing Values Imputation', 'Data', dataField);
    returnDialog('Missing Values Imputation')!.input('Data').input.dispatchEvent(new MouseEvent('click'));
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length == 2, 'cannot find Data columns dialog', 1000);
    okButtonInSelectColumn = returnDialog('Select columns...')?.root
      .getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButtonInSelectColumn.click();
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length == 1, 'error while selectiong data columns', 1000);
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => grok.shell.t.col('AGE')?.stats.missingValueCount === 0,
      'column AGE still contains empty rows!', 9000);
    if (demog.col('AGE')?.stats.missingValueCount !== 0)
      throw new Error('column AGE contains empty rows!');
    if (demog.col('WEIGHT')?.stats.missingValueCount !== 0)
      throw new Error('column WEIGHT contains empty rows!');
    if (demog.col('HEIGHT')?.stats.missingValueCount !== 0)
      throw new Error('column HEIGHT contains empty rows!');
  });

  test('randomData', async () => {
    demog = grok.data.demo.demog();
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 1', 1000);
    isDialogPresent('Random Data');
    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => demog.col('normal') !== null, '"normal" col not added to table', 5000);
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 2', 1000);
    isDialogPresent('Random Data');
    setDialogInputValue('Random Data', 'Distribution', 'log-normal');
    await delay(100);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => demog.col('log-normal') !== null, '"log-normal" col not added to table', 5000);
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 3', 1000);
    setDialogInputValue('Random Data', 'Distribution', 'binomial');
    await delay(100);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => demog.col('binomial') !== null, '"binomial" col not added to table', 5000);
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 4', 1000);
    setDialogInputValue('Random Data', 'Distribution', 'poisson');
    await delay(100);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => demog.col('poisson') !== null, '"poisson" col not added to table', 5000);
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 5', 1000);
    setDialogInputValue('Random Data', 'Distribution', 'uniform');
    await delay(100);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => demog.col('uniform') !== null, '"uniform" col not added to table', 5000);
    isColumnPresent(demog.columns, 'normal');
    isColumnPresent(demog.columns, 'binomial');
    isColumnPresent(demog.columns, 'poisson');
    isColumnPresent(demog.columns, 'uniform');
    isColumnPresent(demog.columns, 'log-normal');
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 6', 1000);
    isDialogPresent('Random Data');
    returnDialog('Random Data')!.input('Show histogram').input.click();
    await awaitCheck(() => Array.from(v.viewers).find((el) => el.type === 'Histogram') !== undefined,
      'Histogram not added', 3500);
    const cancelButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'CANCEL') as HTMLElement;
    cancelButton.click();
    await awaitCheck(() => !checkDialog('Random Data'));
    Array.from(v.viewers).forEach((viewer) => {
      if (viewer.type === 'Histogram')
        throw new Error('Histogram did not disappear after clicking on the "Cancel" button');
    });
    grok.shell.topMenu.find('Ml').find('Tools').find('random Data...').click();
    await awaitCheck(() => checkDialog('Random Data'), 'Dialog is not open 7', 1000);
    isDialogPresent('Random Data');
    returnDialog('Random Data')!.input('Show histogram').input.click();
    await awaitCheck(() => Array.from(v.viewers).find((el) => el.type == 'Histogram') !== undefined,
      'Histogram not added', 3500);
    isViewerPresent(Array.from(v.viewers), 'Histogram');
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton!.click();
    await awaitCheck(() => !checkDialog('Random Data'));
    isViewerPresent(Array.from(v.viewers), 'Histogram');
  });

  test('multivariateAnalysis', async () => {
    demog = df.clone();
    demog.columns.remove('subj');
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => document.querySelector('canvas') !== null, 'cannot load table', 3000);
    grok.shell.topMenu.find('Tools').find('Data Science').find('Multivariate Analysis (PLS)...').click();
    await awaitCheck(() => checkDialog('Multivariate Analysis (PLS)'), 'Dialog is not open', 1000);
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent == 'OK') as HTMLElement;
    if (okButton !== undefined) {
      if (okButton.classList[2] !== 'disabled')
        throw new Error('OK Button does not have disabled class');
    }
    returnDialog('Multivariate Analysis (PLS)')!.input('Features').input.dispatchEvent(new MouseEvent('click'));
    await awaitCheck(() => DG.Dialog.getOpenDialogs().length === 2, 'Select columns dialog did not open', 1000);
    const allButtonInSelectColumns = Array.from(document.querySelectorAll('.d4-link-label'))
      .find((el) => el.textContent === 'All') as HTMLElement;
    const okButtonInSelectColumns = Array.from(returnDialog('Select columns...')!.root
      .querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    if (allButtonInSelectColumns !== undefined) {
      allButtonInSelectColumns.click();
      await delay(100);
    }
    if (okButtonInSelectColumns !== undefined) {
      okButtonInSelectColumns.click();
      await delay(100);
    }
    okButton!.click();
    await awaitCheck(() => Array.from(v.dockManager.rootNode.children).length > 1,
      'Multivariate Analysis is not done', 12000);
  });

  after(async () => {
    grok.shell.windows.simpleMode = false;
  });
});
