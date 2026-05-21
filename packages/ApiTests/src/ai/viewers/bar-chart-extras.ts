import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectPropAndLook, expectRoundTrip} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:682 (BarChartViewer),
// public/js-api/src/interfaces/d4.ts:786 (IBarChartSettings).
// Second pass on Bar chart that pins behaviour the bar-chart-js-api suite,
// the parametric round-trips in src/grid/viewer-set-property.ts, and the
// Dart-side targeted Bar tests don't already cover: the friendly-key `color`
// alias on the typed factory, split + stack JSON envelope shape, the
// barSortType / barSortOrder pair round-trip, choices introspection on
// barSortType, includeNulls boolean round-trip, and the onClick RowGroupAction
// enum. Public JS API only — no canvas geometry.
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
});
