import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import { before, category, test, expect, settled, waitUntil } from '@datagrok-libraries/test/src/test';
import {
  isColumnPresent, isViewerPresent, isDialogPresent, returnDialog,
  setDialogInputValue, checkDialog, checkViewer
} from './gui-utils';


category('GUI: Dialogs', () => {
  let v: DG.TableView;
  let df: DG.DataFrame;
  let demog: DG.DataFrame;

  before(async () => {
    df = grok.data.demo.demog(1000);
  });

  test('anonymize', async () => {
    demog = df.clone();
    v = grok.shell.addTableView(demog);
    await settled(v.grid);
    grok.shell.topMenu.find('Data').find('Anonymize...').click();
    await waitUntil(() => checkDialog('Anonymize Data'), 'Dialog is not open');
    setDialogInputValue('Anonymize Data', 'Number randomization factor', 1);
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await waitUntil(() => !checkDialog('Anonymize Data'), 'dialog was not closed');
  });

  test('aggregateRows', async () => {
    try {
      demog = df.clone();
      v = grok.shell.addTableView(demog);
      await settled(v.grid);
      grok.shell.topMenu.find('Data').find('Aggregate Rows...').click();
      await waitUntil(() => document.querySelector('.grok-pivot-root') !== null, 'the aggregation editor did not open');
      const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok.add-workspace-btn'))
        .find((el) => el.textContent === 'ADD') as HTMLElement;
      okButton.click();
      await waitUntil(() => grok.shell.v.name.toLowerCase() === `${df.name} aggregation`, 'Aggregation table was not created');
      grok.shell.v.close();
    } finally {
      grok.shell.windows.showColumns = false;
    }
  });

  test('cluster', async () => {
    grok.shell.windows.simpleMode = true;
    demog = df.clone();
    v = grok.shell.addTableView(demog);
    await settled(v.grid);
    v.ribbonMenu.find('ML').find('Cluster').find('Cluster Data...').click();
    await waitUntil(() => checkDialog('Cluster Data'), 'Dialog is not open 1');
    isDialogPresent('Cluster Data');
    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await waitUntil(() => demog.col('clusters') !== null, 'cannot find clusters column');
    isColumnPresent(demog.columns, 'clusters');
    v.ribbonMenu.find('ML').find('Cluster').find('Cluster Data...').click();
    await waitUntil(() => checkDialog('Cluster Data'), 'Dialog is not open 2');
    isDialogPresent('Cluster Data');
    returnDialog('Cluster Data')!.input('Show Scatter Plot').input.click();
    await waitUntil(() => checkViewer(Array.from(v.viewers), 'Scatter plot'), 'Scatter plot did not appear');
    isColumnPresent(demog.columns, 'clusters (2)');
    const cancelButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'CANCEL') as HTMLElement;
    cancelButton.click();
    await waitUntil(() => !checkViewer(Array.from(v.viewers), 'Scatter plot'),
      'Scatter plot still displayed');
    Array.from(v.viewers).forEach((viewer) => {
      if (viewer.type === 'Scatter plot')
        throw new Error('Scatter Plot did not disappear after clicking on the "Cancel" button');
    });
    if (demog.col('clusters (2)') !== null)
      throw new Error('cluster (2) column did not disappear after clicking on the "Cancel" button');
    v.ribbonMenu.find('ML').find('Cluster').find('Cluster Data...').click();
    await waitUntil(() => checkDialog('Cluster Data'), 'Dialog is not open 3');
    isDialogPresent('Cluster Data');
    setDialogInputValue('Cluster Data', 'Normalize', 'Z-scores');
    setDialogInputValue('Cluster Data', 'Clusters', 4);
    setDialogInputValue('Cluster Data', 'Metric', 'Manhattan');
    returnDialog('Cluster Data')!.input('Show Scatter Plot').input.click();
    await waitUntil(() => checkViewer(Array.from(v.viewers), 'Scatter plot'),
      'Scatter plot did not appear');
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton.click();
    await waitUntil(() => demog.col('clusters (2)') !== null, 'cannot find clusters (2) column');
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    isColumnPresent(demog.columns, 'clusters (2)');
  });

  test('pca', async () => {
    demog = df.clone();
    v = grok.shell.addTableView(demog);
    await settled(v.grid);
    // EDA registers ML | Analyze | PCA... on the view ribbon; shell.topMenu has no Analyze there.
    const pca = () => v.ribbonMenu.find('ML')?.find('Analyze')?.find('PCA...');
    if (pca() == null) {
      console.log('pca: ML | Analyze | PCA... unavailable (EDA not deployed on this instance) — skipping');
      return;
    }
    pca()!.click();
    await waitUntil(() => checkDialog('PCA'), 'Dialog is not open 1');
    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    // Features is required: without it the input is invalid and PCA cannot be started
    expect(returnDialog('PCA')!.input('Features').validate(), false, 'an empty Features input should be invalid');
    expect(okButton.classList.contains('disabled'), true, 'OK should be disabled without features');
    okButton.click();
    expect(checkDialog('PCA'), true, 'the dialog should stay open when OK is disabled');
    returnDialog('PCA')!.close();
    await waitUntil(() => !checkDialog('PCA'), 'PCA dialog did not close');
    pca()!.click();
    await waitUntil(() => checkDialog('PCA'), 'Dialog is not open 2');
    const featuresField: DG.Column[] = [];
    featuresField.push(demog.col('age') as DG.Column);
    featuresField.push(demog.col('height') as DG.Column);
    featuresField.push(demog.col('weight') as DG.Column);
    setDialogInputValue('PCA', 'Features', featuresField);
    returnDialog('PCA')!.input('Features').input.dispatchEvent(new MouseEvent('click'));
    await waitUntil(() => DG.Dialog.getOpenDialogs().length === 2, 'cannot find Features columns dialog');
    const okButtonInSelectColumn = returnDialog('Select columns...')?.root
      .getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButtonInSelectColumn.click();
    await waitUntil(() => DG.Dialog.getOpenDialogs().length === 1, 'error while selectiong features columns');
    setDialogInputValue('PCA', 'Components', 3);
    returnDialog('PCA')!.input('Center').input.click();
    returnDialog('PCA')!.input('Center').input.click();
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await waitUntil(() => DG.Dialog.getOpenDialogs().length === 0, 'PCA Timeout');
    isColumnPresent(demog.columns, 'PC1');
    isColumnPresent(demog.columns, 'PC2');
    isColumnPresent(demog.columns, 'PC3');
  });

  test('splitColumn', async () => {
    demog = df.clone();
    v = grok.shell.addTableView(demog);
    await settled(v.grid);
    grok.shell.topMenu.find('Data').find('Split Column...').click();
    await waitUntil(() => checkDialog('Split Column'), 'Dialog is not open');
    setDialogInputValue('Split Column', 'Split by', ' ');
    setDialogInputValue('Split Column', 'Prefix', 'TestCol');
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await waitUntil(() => demog.col('TestCol1') !== null && demog.col('TestCol2') !== null,
      'Split Column does not work');
    isColumnPresent(demog.columns, 'TestCol1');
    isColumnPresent(demog.columns, 'TestCol2');
  });

  test('unpivot', async () => {
    demog = df.clone();
    v = grok.shell.addTableView(demog);
    await settled(v.grid);
    grok.shell.topMenu.find('Data').find('Unpivot...').click();
    await waitUntil(() => checkDialog('Unpivot'), 'Dialog is not open');
    const unpivot = returnDialog('Unpivot')!.root;
    const raceInput = Array.from(unpivot.querySelectorAll('td'))
      .find((el) => el.textContent === 'race')!.parentElement!
      .getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;
    const sexInput = Array.from(unpivot.querySelectorAll('td'))
      .find((el) => el.textContent === 'sex')!.parentElement!
      .getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;
    // the selects report a choice to the dialog only through bubbling input and change events
    for (const [input, value] of [[raceInput, 'Copy'], [sexInput, 'Merge']] as [HTMLInputElement, string][]) {
      input.value = value;
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
    }
    const okButton = Array.from(unpivot.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    await waitUntil(() => okButton.classList.contains('enabled'), 'OK stays disabled after choosing the columns');
    okButton.click();
    await waitUntil(() => grok.shell.v.name.toLowerCase().includes('result'), 'result table is not created');
    expect(grok.shell.t.columns.length, 3);
    expect(grok.shell.t.rowCount, 1000);
    grok.shell.v.close();
  });
}, { owner: 'dkovalyov@datagrok.ai' });
