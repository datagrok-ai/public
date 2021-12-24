import {after, before, category, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

let v: DG.TableView;

category('Shell', () => {
  before(async () => {
    v = grok.shell.addTableView(grok.data.demo.demog());
  });

  test('Shell - AddViewer', async () => {


  });

  after(async () => {
    v.close();
    grok.shell.closeTable(v.table!);
  });

});
