import * as DG from 'datagrok-api/dg';
import {category, expect, settled, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectLook, expectPropAndLook, expectRoundTrip, until, withTableView} from '../helpers';

// BarChart JS API extras: friendly-key aliases, sort pair round-trip, choices, includeNulls, onClick enum.
category('AI: Viewers: BarChart extras', () => {
  const v = (): DG.BarChartViewer => DG.Viewer.barChart(demog(), {value: 'age', split: 'race'}) as DG.BarChartViewer;

  test('color friendly-key alias on factory', async () => {
    const c = DG.Viewer.barChart(demog(), {value: 'age', split: 'race', color: 'height'});
    expect(c instanceof DG.BarChartViewer, true);
    expectPropAndLook(c, {colorColumnName: 'height'});
  });

  test('split and stack appear together in look and props', async () => {
    const c = DG.Viewer.barChart(demog(), {value: 'age', split: 'race', stack: 'sex'});
    expectPropAndLook(c, {splitColumnName: 'race', stackColumnName: 'sex'});
  });

  test('barSortType and barSortOrder round-trip together', async () => {
    const c = v();
    c.setOptions({barSortType: 'by value', barSortOrder: 'asc'});
    expectPropAndLook(c, {barSortType: 'by value', barSortOrder: 'asc'});
  });

  test('getProperties choices introspection on barSortType', async () => {
    expectChoices(v(), 'barSortType', ['by category', 'by value']);
  });

  test('includeNulls boolean round-trip via setOptions', async () => {
    const c = v();
    expectRoundTrip(c, {includeNulls: false});
    expectRoundTrip(c, {includeNulls: true});
  });

  test('onClick RowGroupAction enum round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {onClick: 'Filter'});
    expectRoundTrip(c, {onClick: 'Select'});
  });

  test('splitMap month on a datetime split comes back from a layout (public#2562)', async () => {
    await withTableView(demog(1000), async (tv) => {
      const bc = DG.Viewer.barChart(tv.dataFrame, {splitColumnName: 'started', splitMap: 'month'});
      tv.addViewer(bc);
      await settled(bc);
      const layout = tv.saveLayout();
      bc.close();
      tv.loadLayout(layout);
      let restored: DG.Viewer | undefined;
      await until(() => (restored = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.BAR_CHART)) != null);
      expectLook(restored!, {splitColumnName: 'started', splitMap: 'month'});
    });
  });
});
