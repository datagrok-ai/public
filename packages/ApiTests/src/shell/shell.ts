import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const demog = grok.data.demo.demog();

category('Shell', () => {
  test('AddView', async () => {
    const v = grok.shell.addTableView(demog);
    expect(grok.shell.v, v);
    expect(grok.shell.v instanceof DG.TableView, true);
    expect((grok.shell.v as DG.TableView).dataFrame, demog);
    v.close();
    expect(grok.shell.v != v, true);
    expect(grok.shell.table(demog.name), demog);
    grok.shell.closeTable(v.dataFrame);
  }, {owner: 'aparamonov@datagrok.ai'});
});
