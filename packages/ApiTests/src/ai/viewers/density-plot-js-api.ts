import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectLook, expectRoundTrip, look} from '../helpers';

// DensityPlot look round-trips: bins, linearColorScheme, invertColorScheme, colorTransformType, axes, x/y bounds.
category('AI: Viewers: DensityPlot JS API', () => {
  const v = (): DG.Viewer => DG.Viewer.densityPlot(demog(), {x: 'age', y: 'height'});

  test('factory returns typed Viewer with DENSITY_PLOT type', async () => {
    const c = v();
    expect(c.type, DG.VIEWER.DENSITY_PLOT);
    expectLook(c, {xColumnName: 'age', yColumnName: 'height'});
  });

  test('bins integer round-trip within 1-200 range', async () => {
    const c = v();
    for (const n of [120, 1, 200])
      expectRoundTrip(c, {bins: n});
  });

  test('linearColorScheme number-array round-trip', async () => {
    const c = v();
    const scheme = [0xFF112233, 0xFF445566, 0xFFAABBCC];
    c.setOptions({linearColorScheme: scheme});
    expectArray(look(c)['linearColorScheme'] as number[], scheme);
  });

  test('showXAxis/showYAxis/showColorScale combined bool round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {showXAxis: false, showYAxis: false, showColorScale: false});
    expectRoundTrip(c, {showXAxis: true, showYAxis: true, showColorScale: true});
  });

  test('x/yMin x/yMax numeric round-trip', async () => {
    expectRoundTrip(v(), {xMin: 10, xMax: 80, yMin: 120, yMax: 200});
  });

  test('invertColorScheme + colorTransformType logarithmic round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {invertColorScheme: true, colorTransformType: 'logarithmic'});
    expectRoundTrip(c, {invertColorScheme: false, colorTransformType: 'linear'});
  });

  test('getProperties choices introspection on colorTransformType', async () => {
    expectChoices(v(), 'colorTransformType', ['linear', 'logarithmic']);
  });
}, {owner: 'agolovko@datagrok.ai'});
