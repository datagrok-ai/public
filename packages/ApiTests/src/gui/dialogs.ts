import {after, before, awaitCheck, category, test, delay, expect} from '@datagrok-libraries/utils/src/test';
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
    await awaitCheck((() => {return checkViewer(Array.from(v.viewers), 'Scatter plot')}), 'Scatter plot did not appear', 2500); 
    isColumnPresent(demog.columns, 'clusters (2)');

    const cancelButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'CANCEL') as HTMLElement;
    cancelButton.click();

    await awaitCheck((() => {return !checkViewer(Array.from(v.viewers), 'Scatter plot')}), 'Scatter plot still displayed', 1000); 

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
    await awaitCheck((() => {return checkViewer(Array.from(v.viewers), 'Scatter plot')}), 'Scatter plot did not appear', 2500); 

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
    await awaitCheck((() => {return grok.shell.t.col('AGE')?.stats.missingValueCount == 0}), 'column AGE still contains empty rows!', 9000);

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
    await awaitCheck((() => {return DG.Dialog.getOpenDialogs().length == 0}), 'PCA Timeout', 2500);

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
    await awaitCheck((() => {return grok.shell.t.col('normal') != null}), '"normal" col not added to table', 2500);

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    isDialogPresent('Random Data');
    setDialogInputValue('Random Data', 'Distribution', 'log-normal'); await delay(200);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click(); 
    await awaitCheck((() => {return grok.shell.t.col('log-normal') != null}), '"log-normal" col not added to table', 2500);

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    setDialogInputValue('Random Data', 'Distribution', 'binomial'); await delay(200);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck((() => {return grok.shell.t.col('binomial') != null}), '"binomial" col not added to table', 2500);

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    setDialogInputValue('Random Data', 'Distribution', 'poisson'); await delay(200);
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); 
    await awaitCheck((() => {return grok.shell.t.col('poisson') != null}), '"poisson" col not added to table',2500);

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    setDialogInputValue('Random Data', 'Distribution', 'uniform'); await delay(200);
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton!.click(); 
    await awaitCheck((() => {return grok.shell.t.col('uniform') != null}),'"uniform" col not added to table', 2500);

    isColumnPresent(demog.columns, 'normal');
    isColumnPresent(demog.columns, 'binomial');
    isColumnPresent(demog.columns, 'poisson');
    isColumnPresent(demog.columns, 'uniform');
    isColumnPresent(demog.columns, 'log-normal');

    grok.shell.topMenu.find('Tools').find('Data Science').find('Random Data...').click(); 
    await awaitCheck(() => {return checkDialog('Random Data');});
    isDialogPresent('Random Data');
    returnDialog('Random Data')!.input('Show histogram').input.click();
    await awaitCheck((() => {return Array.from(v.viewers).find((el) => el.type == 'Histogram') != undefined}), 'Histogram not added', 3500);

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
    returnDialog('Random Data')!.input('Show histogram').input.click();
    await awaitCheck((() => {return Array.from(v.viewers).find((el) => el.type == 'Histogram') != undefined}), 'Histogram not added', 3500);
    isViewerPresent(Array.from(v.viewers), 'Histogram');
    okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton!.click(); 
    await awaitCheck(() => {return !checkDialog('Random Data');});

    isViewerPresent(Array.from(v.viewers), 'Histogram');

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });  

  test('dialogs.multivariateAnalysis', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Multivariate Analysis (PLS)...').click();
    await awaitCheck((() => {return document.querySelector('.d4-dialog') != undefined}), 'Dialog did not open', 1000);

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent == 'OK') as HTMLElement;
    
    if (okButton != undefined) {
      if (okButton.classList[2] != 'disabled')
        throw 'OK Button does not have disabled class';
    }

    returnDialog('Multivariate Analysis (PLS)')!.input('Features').input.dispatchEvent(new MouseEvent('click')); 
    await awaitCheck((() => {return DG.Dialog.getOpenDialogs().length == 2}), 'Select columns dialog did not open', 1000);

    let allButtonInSelectColumns = Array.from(document.querySelectorAll('.d4-link-label'))
    .find((el) => el.textContent === 'All') as HTMLElement;

    let okButtonInSelectColumns = Array.from(returnDialog('Select columns...')!.root.querySelectorAll('.ui-btn.ui-btn-ok'))
    .find((el) => el.textContent === 'OK') as HTMLElement;

    if (allButtonInSelectColumns != undefined) {
      allButtonInSelectColumns.click(); await delay(200);
    }
    if (okButtonInSelectColumns != undefined) {
      okButtonInSelectColumns.click(); await delay(200);
    }
    okButton!.click(); 
    await awaitCheck((() => {return Array.from(v.dockManager.rootNode.children).length > 1}), 'Multivariate Analysis is not done', 12000);

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });
  
  test('dialogs.spe', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Tools').find('Data Science').find('Stochastic Proximity Embedding...').click();
    await awaitCheck(() => {return checkDialog('SPE');});

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    //WIP, SPE Feature currently not working (GROK-10250)

    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  }, {skipReason: 'GROK-10250'});  

  test('dialogs.splitColumn', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});

    grok.shell.topMenu.find('Data').find('Split Column...').click();
    await awaitCheck(() => {return checkDialog('Text to columns');});

    setDialogInputValue('Text to columns', 'Split by', ' '); await delay(200);
    setDialogInputValue('Text to columns', 'Prefix', 'TestCol'); await delay(200);

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck((() => demog.col('TestCol1') != undefined && demog.col('TestCol2') != undefined), '', 3000);

    isColumnPresent(demog.columns, 'TestCol1');
    isColumnPresent(demog.columns, 'TestCol2');
    
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });
  
  test('dialogs.unpivot', async () => {
    let v: DG.TableView;
    const demog = grok.data.demo.demog(1000);
    v = grok.shell.addTableView(demog);
    await awaitCheck(() => {return grok.shell.v == v});
    
    grok.shell.topMenu.find('Data').find('Unpivot...').click(); 
    await awaitCheck(() => {return checkDialog('Unpivot');});

    let raceInput = Array.from(document.querySelectorAll('td')).find((el) => el.textContent == 'race')!.parentElement!.getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;
    let sexInput = Array.from(document.querySelectorAll('td')).find((el) => el.textContent == 'sex')!.parentElement!.getElementsByClassName('ui-input-editor')[0] as HTMLInputElement;

    raceInput.value = 'Copy';
    sexInput.value = 'Merge';

    let okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
        .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck((() => {return grok.shell.v.name == 'result'}), 'result table is not created', 2000)

    expect(grok.shell.t.columns.length, 3);
    expect(grok.shell.t.rowCount, 1000);
    
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });  
});
