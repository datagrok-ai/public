import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let demog = grok.data.demo.demog();

category('Shell', () => {
  before(async () => {
    grok.shell.closeAll();
  });

  test('AddView', async () => {
    let v = grok.shell.addTableView(demog);
    expect(grok.shell.v, v);
    expect(grok.shell.v instanceof DG.TableView, true);
    expect((grok.shell.v as DG.TableView).dataFrame, demog);
    v.close();
    expect(grok.shell.v, null);
    expect(grok.shell.table(demog.name), demog);
    grok.shell.closeTable(v.dataFrame);
  });

  after(async () => {
 //   v.close();
   // grok.shell.closeTable(v.table!);
  });

});
