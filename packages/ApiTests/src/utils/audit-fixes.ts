import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

// Regressions for the defects fixed after the 2026-09 JS API audit (core/docs/reviews/js-api-audit/defects.md).
category('JS API: audit fixes', () => {
  test('DataFrame.onDataChanged fires on column and row changes', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.TYPE.INT, 'x', [1, 2, 3])]);
    let fired = 0;
    const sub = df.onDataChanged.subscribe(() => fired++);
    df.columns.addNewInt('y');
    df.rows.addNew([4, 5]);
    sub.unsubscribe();
    expect(fired >= 2, true, `expected column-add and row-add events, got ${fired}`);
  });

  test('Rect.getGridPart divides the height by yCount', async () => {
    const part = new DG.Rect(0, 0, 100, 40).getGridPart(2, 2, 1, 1);
    expect(part.x, 50);
    expect(part.y, 20);
    expect(part.width, 50);
    expect(part.height, 20);
  });

  test('Rect.fromCenterSize is static', async () => {
    const r = DG.Rect.fromCenterSize(50, 50, 20, 10);
    expect(r.x, 40);
    expect(r.y, 45);
    expect(r.width, 20);
    expect(r.height, 10);
  });

  test('FormulaLinesHelper.removeAt removes only the requested items', async () => {
    const lines = grok.data.demo.demog(10).meta.formulaLines;
    lines.clear();
    lines.addAll([{formula: '${age} = 1'}, {formula: '${age} = 2'}, {formula: '${age} = 3'}]);
    lines.removeAt(0);
    expect(lines.items.length, 2);
    expect(lines.items[0].formula, '${age} = 2');
    lines.removeAt(0, 2);
    expect(lines.items.length, 0);
  });

  test('Color.hexToPercentRgb scales by 255', async () => {
    expectArray(DG.Color.hexToPercentRgb('#ffffff')!, [1, 1, 1, 0.3]);
  });

  test('ui.iconImage resolves a bare file name under /images', async () => {
    const icon = ui.iconImage('logo', 'logo.png');
    expect(icon.style.backgroundImage.includes('/images/logo.png'), true, icon.style.backgroundImage);
  });

  test('viewer.props reports only existing properties for "in"', async () => {
    const viewer = DG.Viewer.scatterPlot(grok.data.demo.demog(10));
    try {
      expect('xColumnName' in viewer.props, true);
      expect('noSuchProperty' in viewer.props, false);
    } finally {
      viewer.detach();
    }
  });

  test('TableQueryBuilder join helpers keep the alias', async () => {
    const query = DG.TableQueryBuilder.from('orders').leftJoin('customers', ['customer_id'], ['id'], 'c').build();
    expect(query.joins.length, 1);
    expect(query.joins[0].rightTableAlias, 'c');
    expect(query.joins[0].joinType, DG.JOIN_TYPE.LEFT);
  });

  test('GroupByBuilder.getGroups returns a plain object', async () => {
    const groups = grok.data.demo.demog(100).groupBy(['sex']).getGroups();
    expect(Object.keys(groups).length > 0, true);
    expect(Object.values(groups)[0] instanceof DG.DataFrame, true);
  });
});
