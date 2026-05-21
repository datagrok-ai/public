import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectCleared, expectRoundTrip, findProp, look} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:219 (DG.Viewer.histogram),
// public/js-api/src/viewer.ts:671 (HistogramViewer),
// public/js-api/src/interfaces/d4.ts:1733 (IHistogramSettings),
// public/js-api/src/entities/property.ts:121 (Property — choices,
// columnTypeFilter).
// Histogram-only JSON-shape coverage that the Dart parametric harness can't
// reach through its boolean/midpoint/column iterator strategy: clearing a
// *ColumnName via empty string, formula filter strings, array-typed look
// values (linearColorScheme), and getProperties() introspection on
// columnTypeFilter and colorAggrType choices. All assertions read state via
// getOptions(true).look — no first-paint geometry.
category('AI: Viewers: Histogram extras', () => {
  const v = (): DG.HistogramViewer => DG.Viewer.histogram(demog(), {value: 'age'});

  test('clear splitColumnName via setOptions empty string', async () => {
    const c = v();
    expectRoundTrip(c, {splitColumnName: 'race'});
    c.setOptions({splitColumnName: ''});
    expectCleared(look(c)['splitColumnName']);
    expectCleared(c.props['splitColumnName']);
    expect(look(c)['splitColumnName'], c.props['splitColumnName']);
  });

  test('filter formula round-trip', async () => {
    expectRoundTrip(v(), {filter: '${age} > 20'});
  });

  test('linearColorScheme number-array round-trip', async () => {
    const c = v();
    const scheme = [0xFF112233, 0xFFAABBCC];
    c.setOptions({linearColorScheme: scheme});
    expectArray(look(c)['linearColorScheme'] as number[], scheme);
  });

  test('getProperties columnTypeFilter introspection', async () => {
    const c = v();
    expect(findProp(c, 'valueColumnName')!.columnTypeFilter as string, 'numerical');
    expect(findProp(c, 'splitColumnName')!.columnTypeFilter as string, 'categorical');
    expect(findProp(c, 'colorColumnName')!.columnTypeFilter as string, 'numerical');
  });

  test('getProperties choices introspection on colorAggrType', async () => {
    const aggrProp = findProp(v(), 'colorAggrType');
    expect(aggrProp != null, true);
    expect(Array.isArray(aggrProp!.choices), true);
    expect(aggrProp!.choices.length > 0, true);
  });

  test('aggTooltipColumns string round-trip', async () => {
    expectRoundTrip(v(), {aggTooltipColumns: 'avg(age), count(race)'});
  });
});
