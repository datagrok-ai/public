import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getOptions} from './utils';


category('Word cloud', () => {
  const TYPE = 'Word cloud';
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog(100);
    tv = grok.shell.addTableView(df);
  });

  test('Legacy columnColumnName layout', async () => {
    const viewer = tv.addViewer(TYPE, {wordColumnName: 'race'});
    await awaitCheck(() => viewer.getOptions(true).look['wordColumnName'] === 'race', 'option not applied', 3000);
    const json = tv.saveLayout().toJson().replace(/"wordColumnName"/g, '"columnColumnName"');
    viewer.close();

    tv.loadLayout(DG.ViewLayout.fromJson(json));
    let restored: DG.Viewer | undefined;
    await awaitCheck(() => (restored = Array.from(tv.viewers).find((v) => v.type === TYPE)) != null,
      'layout did not restore the viewer', 3000);
    expect((await getOptions(restored!)).wordColumnName, 'race');
  });

  test('Long words', async () => {
    const words = Array.from({length: 5}, (_, i) => `w${i} ${'x'.repeat(5000)}`);
    const longDf = DG.DataFrame.fromColumns([DG.Column.fromStrings('text', words.concat(words))]);
    const longTv = grok.shell.addTableView(longDf);
    try {
      const viewer = longTv.addViewer(TYPE, {wordColumnName: 'text'});
      const status = () => DG.Widget.find(viewer.root)!.getWidgetStatus();
      await awaitCheck(() => status().values?.['words'] === words.length, 'long words were not laid out', 5000);
      expect(status().hitAreas[`word "${words[0]}"`] !== undefined, true);
    }
    finally {
      longTv.close();
      grok.shell.closeTable(longDf);
    }
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
