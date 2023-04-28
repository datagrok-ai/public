import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import wu from 'wu';


category('View: Docking', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  test('dock', async () => {
    expect(wu(tv.viewers).find((v) => v.tags['test']) == undefined);
    const viewer = df.plot.scatter();
    viewer.tags['test'] = 'true';
    tv.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN);
    expect(wu(tv.viewers).find((v) => v.tags['test'])?.tags['test'], 'true');
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
