import {after, awaitCheck, before, category, delay, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isColumnPresent, isDialogPresent, setDialogInputValue} from './gui-utils';


category('GUI', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('dialogs.addNewColumn', async () => {
    const addNewColumnBtn = document.getElementsByClassName('svg-add-new-column')[0] as HTMLElement;
    addNewColumnBtn.click(); 
    await awaitCheck(() => {return document.querySelector('.d4-dialog') != undefined;});
    setDialogInputValue('Add New Column', 'Name', 'TestNewColumn');
    setDialogInputValue('Add New Column', 'Expression', '${AGE}+5');
    await delay(500);
    let okButton = document.getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
    okButton!.click(); 
    await awaitCheck(() => {return demog.col('TestNewColumn') != undefined;});
    isColumnPresent(demog.columns, 'TestNewColumn');
    addNewColumnBtn.click(); 
    await awaitCheck(() => {return document.querySelector('.d4-dialog') != undefined;});
    isDialogPresent('Add New Column');
    setDialogInputValue('Add New Column', 'Name', 'TestNewStringColumn');
    setDialogInputValue('Add New Column', 'Expression', 'TestValue');
    await delay(500);
    okButton = document.getElementsByClassName('ui-btn ui-btn-ok')[0] as HTMLElement;
    okButton!.click(); 
    await awaitCheck(() => {return demog.col('TestNewStringColumn') != undefined;});
    isColumnPresent(demog.columns, 'TestNewStringColumn');
  });

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });
});
