import {after, before, awaitCheck, category, test, isDialogPresent} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {setDialogInputValue} from './gui-utils';


category('GUI', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('dialogs.anonymize', async () => {
    grok.shell.topMenu.find('Data').find('Anonymize...').click(); 
    await awaitCheck(() => {return isDialogPresent('Anonymize Data');});
    setDialogInputValue('Anonymize Data', 'Number randomization factor', 1);
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    await awaitCheck(() => {return !isDialogPresent('Anonymize Data');});
  });

  after(async () => {
    v.close();
  });
});
