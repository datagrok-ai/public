import {after, before, category, delay, awaitCheck, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue, waitForElement, checkDialog, checkViewer} from './gui-utils';

category('GUI: Cluster', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('dialogs.cluster', async () => {
    grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').click(); 
    await awaitCheck(() => {return checkDialog('Cluster Data')}); 
    
    isDialogPresent('Cluster Data');

    let okButton:HTMLElement | undefined;
    let button;
    for (let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++) {
      button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
      if (button.innerText == 'OK')
        okButton = button;
    }
    okButton!.click();
    await awaitCheck(() => {return demog.col('clusters') != undefined}); 
    isColumnPresent(demog.columns, 'clusters');

    grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').click(); 
    await awaitCheck(() => {return checkDialog('Cluster Data')}); 
    isDialogPresent('Cluster Data');

    returnDialog('Cluster Data')!.input('Show scatter plot').input.click(); await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return checkViewer(Array.from(v.viewers), 'Scatter plot')}); 
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    isColumnPresent(demog.columns, 'clusters (2)');

    let cancelButton:HTMLElement;
    for (let i=0; i<document.getElementsByClassName('ui-btn ui-btn-ok').length; i++) {
      button = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
      if (button.innerText == 'CANCEL') {
        cancelButton = document.getElementsByClassName('ui-btn ui-btn-ok')[i] as HTMLElement;
        cancelButton.click();
      }
    }

    await awaitCheck(() => {return !checkViewer(Array.from(v.viewers), 'Scatter plot')}); 

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
      if (Array.from(v.viewers)[i].type == 'Scatter plot')
        throw 'Scatter Plot did not disappear after clicking on the "Cancel" button';
    }

    if (demog.columns.byName('clusters (2)') != null)
      throw 'cluster (2) column did not disappear after clicking on the "Cancel" button';

    grok.shell.topMenu.find('Tools').find('Data Science').find('Cluster...').click(); 
    await awaitCheck(() => {return checkDialog('Cluster Data')}); 
    isDialogPresent('Cluster Data');

    setDialogInputValue('Cluster Data', 'Normalize', 'Z-scores'); await delay(100);
    setDialogInputValue('Cluster Data', 'Clusters', 4); await delay(100);
    setDialogInputValue('Cluster Data', 'Metric', 'Manhattan'); await delay(100);
    returnDialog('Cluster Data')!.input('Show scatter plot').input.click(); await delay(2000); // not enough timeout inside awaitCheck(). I will remove it after the parameterization of the function
    await awaitCheck(() => {return checkViewer(Array.from(v.viewers), 'Scatter plot')});

    okButton = document.getElementsByClassName('ui-btn ui-btn-ok enabled')[0] as HTMLElement;
    okButton!.click(); 
    await awaitCheck(() => {return demog.col('clusters (2)') != undefined}); 
    
    isViewerPresent(Array.from(v.viewers), 'Scatter plot');
    isColumnPresent(demog.columns, 'clusters (2)');
  });

  after(async () => {
    grok.shell.closeAll();
  });
});
