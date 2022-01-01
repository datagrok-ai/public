import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { after, before, category, expect, test } from '../test';

category('Layouts', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog();
    tv = grok.shell.addTableView(df);
  });

  test('TableView.saveLayout()', async () => {
    expect(tv.saveLayout() instanceof DG.ViewLayout, true);
  });

  test('ViewLayout.toJson()', async () => {
    const layout = tv.saveLayout();
    const json = layout.toJson();
    expect(typeof json, 'string');
    expect(JSON.parse(json).id, layout.id);
  });

  test('ViewLayout.fromJson(string)', async () => {
    const layout = tv.saveLayout();
    const json = layout.toJson();
    expect(DG.ViewLayout.fromJson(json) instanceof DG.ViewLayout, true);
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
