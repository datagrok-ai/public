import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';


category('SankeyViewer', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;
  let viewer: DG.Viewer;

  before(async () => {
    df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('source', ['A', 'B', 'A', 'B']),
      DG.Column.fromStrings('target', ['B', 'C', 'C', 'C']),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'value', [5, 3, 0, null]),
    ]);
    tv = grok.shell.addTableView(df);
    viewer = tv.addViewer('Sankey');
    await awaitCheck(() => viewer.root.querySelector('svg') !== null, 'Sankey has not rendered', 3000);
  });

  test('Filtering out all rows', async () => {
    df.filter.setAll(false);
    await awaitCheck(() => viewer.root.querySelector('svg') === null, 'Sankey still draws filtered-out rows', 3000);
    df.filter.setAll(true);
    await awaitCheck(() => viewer.root.querySelector('svg') !== null, 'Sankey has not rendered', 3000);
  });

  test('Filtering to rows without weight', async () => {
    df.filter.init((i) => i >= 2);
    await awaitCheck(() => viewer.root.querySelector('svg') === null, 'Sankey draws links without weight', 3000);
    expect(viewer.root.querySelectorAll('rect[height="NaN"]').length, 0);
    df.filter.setAll(true);
    await awaitCheck(() => viewer.root.querySelector('svg') !== null, 'Sankey has not rendered', 3000);
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
