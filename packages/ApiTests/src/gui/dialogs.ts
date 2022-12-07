import {after, before, awaitCheck, category, test, delay} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog,
    setDialogInputValue, checkDialog, checkViewer} from './gui-utils';


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
});
