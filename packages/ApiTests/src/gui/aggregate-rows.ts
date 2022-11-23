import {after, before, awaitCheck, category, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


category('GUI', () => {
  let v: DG.TableView;
  const demog = grok.data.demo.demog(1000);

  before(async () => {
    v = grok.shell.addTableView(demog);
  });

  test('dialogs.aggregateRows', async () => {
    grok.shell.topMenu.find('Data').find('Aggregate Rows...').click(); 
    await awaitCheck(() => {return document.querySelector('.grok-pivot') != undefined;}); 
    const okButton = Array.from(document.querySelectorAll('.ui-btn.ui-btn-ok'))
      .find((el) => el.textContent === 'OK') as HTMLElement;
    okButton.click();
    try {
      await awaitCheck(() => {return grok.shell.v.name == 'result';});
      grok.shell.v.close();
    } catch (e) {
      throw new Error('Aggregation table was not created');
    } finally {
      grok.shell.windows.showColumns = false;
    }
  });

  after(async () => {
    v.close();
    grok.shell.tables.forEach((t) => grok.shell.closeTable(t));
  });
});
