import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectPropAndLook, expectRoundTrip} from '../helpers';

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
});
