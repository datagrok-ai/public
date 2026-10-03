import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';


category('Sunburst', () => {
  let df: DG.DataFrame;
  let tv: DG.TableView;

  before(async () => {
    df = grok.data.demo.demog(100);
    tv = grok.shell.addTableView(df);
  });

  test('Row reordering in another grid', async () => {
    const sunburst = tv.addViewer('Sunburst');
    await awaitCheck(() => (sunburst.props.hierarchyColumnNames?.length ?? 0) > 0, 'hierarchy not set', 5000);
    grok.shell.o = sunburst;
    const hierarchy = [...sunburst.props.hierarchyColumnNames];

    const names = DG.DataFrame.fromColumns([DG.Column.fromStrings('name', ['a', 'b', 'c'])]);
    const grid = DG.Viewer.grid(names, {allowRowReordering: true});
    grid.root.style.width = '200px';
    grid.root.style.height = '200px';
    tv.dockManager.dock(grid.root, 'right');
    await delay(500);

    const errors: string[] = [];
    const onError = (e: ErrorEvent) => errors.push(e.message);
    window.addEventListener('error', onError);
    try {
      const bounds = grid.cell('name', 0).bounds;
      const rect = grid.overlay.getBoundingClientRect();
      const x = rect.left + bounds.x + bounds.width / 2;
      const y = rect.top + bounds.y + bounds.height / 2;
      const mouse = (dy: number) =>
        ({bubbles: true, cancelable: true, clientX: x, clientY: y + dy, button: 0, buttons: 1, view: window});
      grid.overlay.dispatchEvent(new MouseEvent('mousedown', mouse(0)));
      document.dispatchEvent(new MouseEvent('mousemove', mouse(10)));
      document.dispatchEvent(new MouseEvent('mousemove', mouse(2 * bounds.height)));
      document.dispatchEvent(new MouseEvent('mouseup', mouse(2 * bounds.height)));
      await delay(500);
    } finally {
      window.removeEventListener('error', onError);
    }

    expect(errors.length, 0, errors.join('; '));
    expect(JSON.stringify(sunburst.props.hierarchyColumnNames), JSON.stringify(hierarchy));
  });

  after(async () => {
    tv.close();
    grok.shell.closeTable(df);
  });
});
