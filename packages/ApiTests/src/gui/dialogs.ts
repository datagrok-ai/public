import {after, before, awaitCheck, category, test, delay} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, checkDialog, checkViewer, isErrorBallon} from './gui-utils';


category('GUI: Dialogs', () => {

  test('dialogs.anonymize', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Data').find('Anonymize...').click(); 
    await awaitCheck(() => {return checkDialog('Anonymize Data');});
    setDialogInputValue('Anonymize Data', 'Number randomization factor', 1);
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => {return !checkDialog('Anonymize Data');});    
  });

  test('dialogs.aggregateRows', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Data').find('Aggregate Rows...').click(); 
    await awaitCheck(() => {return document.querySelector('.grok-pivot') != undefined;}); 
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    try {
      await awaitCheck(() => {return grok.shell.v.name == 'result';});
    } catch (e) {
      throw new Error('Aggregation table was not created');
    } finally {
      grok.shell.windows.showColumns = false;
    }
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });

  test('dialogs.cluster', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').click(); 
    await awaitCheck(() => {return checkDialog('Cluster Data');}); 
    isDialogPresent('Cluster Data');
    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
   
    await awaitCheck(() => {return demog.col('clusters') != undefined;});    
    isColumnPresent(demog.columns, 'clusters');

    grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').click(); 
    await awaitCheck(() => {return checkDialog('Cluster Data');}); 
    isDialogPresent('Cluster Data');

    returnDialog('Cluster Data')!.input('Show scatter plot').input.click();    
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return checkViewer(Array.from(v.viewers), 'Scatter plot');}); 
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    isColumnPresent(demog.columns, 'clusters (2)');

    const cancelButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'CANCEL') as HTMLElement;
    cancelButton.click();

    await awaitCheck(() => {return !checkViewer(Array.from(v.viewers), 'Scatter plot');}); 

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
      if (Array.from(v.viewers)[i].type == 'Scatter plot')
        throw new Error('Scatter Plot did not disappear after clicking on the "Cancel" button');
    }

    if (demog.columns.byName('clusters (2)') != null)
      throw new Error('cluster (2) column did not disappear after clicking on the "Cancel" button');

    grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').click(); 
    await awaitCheck(() => {return checkDialog('Cluster Data');}); 
    isDialogPresent('Cluster Data');

    setDialogInputValue('Cluster Data', 'Normalize', 'Z-scores'); await delay(100);
    setDialogInputValue('Cluster Data', 'Clusters', 4); await delay(100);
    setDialogInputValue('Cluster Data', 'Metric', 'Manhattan'); await delay(100);
    returnDialog('Cluster Data')!.input('Show scatter plot').input.click();    
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return checkViewer(Array.from(v.viewers), 'Scatter plot');});

    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); 
    await awaitCheck(() => {return demog.col('clusters (2)') != undefined;}); 
    
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    isColumnPresent(demog.columns, 'clusters (2)');

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });

  test('dialogs.missingValuesImputation', async () => {
    let v:DG.TableView;
    const demog = grok.data.loadTable('https://dev.datagrok.ai/demo/demog.csv');
    v = grok.shell.addTableView(await demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Missing Values Imputation...').click();
    await awaitCheck(() => {return checkDialog('Missing Values Imputation');});
    isDialogPresent('Missing Values Imputation');

    let imputeField:DG.Column[] = [];
    imputeField.push((await demog).columns.byName('AGE'));
    imputeField.push((await demog).columns.byName('HEIGHT'));
    imputeField.push((await demog).columns.byName('WEIGHT'));
    setDialogInputValue('Missing Values Imputation', 'Impute', imputeField);
    returnDialog('Missing Values Imputation')!.input('Impute').input.dispatchEvent(new MouseEvent('click')); 
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 2})
    let okButtonInSelectColumn = returnDialog('Select columns...')?.root.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButtonInSelectColumn.click(); 
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 1})

    let dataField:DG.Column[] = [];
    dataField.push((await demog).columns.byName('SEX'));
    dataField.push((await demog).columns.byName('RACE'));
    dataField.push((await demog).columns.byName('DIS_POP'));
    dataField.push((await demog).columns.byName('DEMOG'));
    setDialogInputValue('Missing Values Imputation', 'Data', dataField);
    returnDialog('Missing Values Imputation')!.input('Data').input.dispatchEvent(new MouseEvent('click')); 
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 2})
    okButtonInSelectColumn = returnDialog('Select columns...')?.root.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButtonInSelectColumn.click(); 
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 1})

    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await delay(8000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function. Impute missed values takes a long time

    if (((await demog).col('AGE')?.stats.missingValueCount) != 0)
      throw 'column AGE contains empty rows!';
    if (((await demog).col('WEIGHT')?.stats.missingValueCount) != 0)
      throw 'column WEIGHT contains empty rows!';
    if (((await demog).col('HEIGHT')?.stats.missingValueCount) != 0)
      throw 'column HEIGHT contains empty rows!';
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });

  test('dialogs.pca', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Principal Component Analysis...').click();
    await awaitCheck(() => {return checkDialog('PCA')});

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => {return !checkDialog('PCA')});

    isErrorBallon('Errors calling PCA: features: Value not defined.');

    grok.shell.topMenu.find('Tools').find('Data Science').find('Principal Component Analysis...').click(); 
    await awaitCheck(() => {return checkDialog('PCA')});

    let featuresField:DG.Column[] = [];
    featuresField.push(demog.columns.byName('age'));
    featuresField.push(demog.columns.byName('height'));
    featuresField.push(demog.columns.byName('weight'));
    setDialogInputValue('PCA', 'Features', featuresField);

    returnDialog('PCA')!.input('Features').input.dispatchEvent(new MouseEvent('click'));
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 2});
    const okButtonInSelectColumn = returnDialog('Select columns...')?.root.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButtonInSelectColumn.click(); 
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 1});

    setDialogInputValue('PCA', 'Components', 3); await delay(200);
    returnDialog('PCA')!.input('Center').input.click(); 
    await awaitCheck(() => {return returnDialog('PCA')!.input('Center').value == false});

    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
    .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton!.click(); 
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function. PCA function takes more time
    await awaitCheck(() => {return DG.Dialog.getOpenDialogs().length == 0});

    isColumnPresent(demog.columns, 'PCA0');
    isColumnPresent(demog.columns, 'PCA1');
    isColumnPresent(demog.columns, 'PCA2');

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });

  test('dialogs.randomData', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click();
    await awaitCheck(() => {return checkDialog('Random Data');});
    isDialogPresent('Random Data');

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return grok.shell.t.col('normal') != null});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    isDialogPresent('Random Data');
    setDialogInputValue('Random Data', 'Distribution', 'log-normal'); await delay(200);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click(); 
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return grok.shell.t.col('log-normal') != null});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    setDialogInputValue('Random Data', 'Distribution', 'binomial'); await delay(200);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return grok.shell.t.col('binomial') != null});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    setDialogInputValue('Random Data', 'Distribution', 'poisson'); await delay(500);
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); 
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return grok.shell.t.col('poisson') != null});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    setDialogInputValue('Random Data', 'Distribution', 'uniform'); await delay(200);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton!.click(); 
    await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return grok.shell.t.col('uniform') != null});

    isColumnPresent(demog.columns, 'normal');
    isColumnPresent(demog.columns, 'binomial');
    isColumnPresent(demog.columns, 'poisson');
    isColumnPresent(demog.columns, 'uniform');
    isColumnPresent(demog.columns, 'log-normal');

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    isDialogPresent('Random Data');
    returnDialog('Random Data')!.input('Show histogram').input.click(); await delay(3000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return Array.from(v.viewers).find((el) => el.type == 'Histogram') != undefined});

    let cancelButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
    .find((el) => el.textContent === 'CANCEL') as HTMLElement;
    cancelButton.click(); 
    await awaitCheck(() => {return !checkDialog('Random Data');});

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
      if (Array.from(v.viewers)[i].type == 'Histogram')
        throw 'Histogram did not disappear after clicking on the "Cancel" button';
    }

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    isDialogPresent('Random Data');
    returnDialog('Random Data')!.input('Show histogram').input.click(); await delay(3000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return Array.from(v.viewers).find((el) => el.type == 'Histogram') != undefined});
    isViewerPresent(Array.from(v.viewers), 'Histogram');
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton!.click(); 
    await awaitCheck(() => {return !checkDialog('Random Data');});

    isViewerPresent(Array.from(v.viewers), 'Histogram');

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });  
});
