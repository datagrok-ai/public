import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { after, before, category, expect, test } from '../test';

category('Layouts', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;
  let baseLayout: DG.ViewLayout;

  before(async () => {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  test('TableView.saveLayout()', async () => {
    baseLayout = tv.saveLayout();
    expect(baseLayout instanceof DG.ViewLayout, true);
  });

  test('ViewLayout.toJson()', async () => {
    const json = baseLayout.toJson();
    expect(typeof json, 'string');
    expect(JSON.parse(json).id, baseLayout.id);
  });

  test('ViewLayout.fromJson(string)', async () => {
    expect(DG.ViewLayout.fromJson(baseLayout.toJson()) instanceof DG.ViewLayout, true);
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
