import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, look} from '../helpers';

// Regression coverage for GROK-19760 (also GROK-19581): Histogram survives invalid valueMin/valueMax.
category('AI: GROK-19760: Histogram invalid valueMin/valueMax', () => {
  function assertLive(v: DG.Viewer, rowCount: number): void {
    expect(v.dataFrame != null, true);
    expect(v.dataFrame!.rowCount, rowCount);
    const l = look(v);
    expect(v.props['valueMin'] == l['valueMin'] || (v.props['valueMin'] == null && l['valueMin'] == null), true);
    expect(v.props['valueMax'] == l['valueMax'] || (v.props['valueMax'] == null && l['valueMax'] == null), true);
  }

  test('inverted range does not throw and survives', async () => {
    const df = demog();
    const v = df.plot.histogram({value: 'age'});
    expectNoThrow(() => v.setOptions({valueMin: 100, valueMax: 0}));
    assertLive(v, df.rowCount);
  });

  test('NaN valueMin does not throw', async () => {
    const df = demog();
    const v = df.plot.histogram({value: 'age'});
    expectNoThrow(() => v.setOptions({valueMin: NaN, valueMax: 50}));
    assertLive(v, df.rowCount);
  });

  test('zero range survives, then a valid range lands cleanly', async () => {
    const df = demog();
    const v = df.plot.histogram({value: 'age'});
    expectNoThrow(() => v.setOptions({valueMin: 50, valueMax: 50}));
    assertLive(v, df.rowCount);
    expectNoThrow(() => v.setOptions({valueMin: 0, valueMax: 100}));
    expectLook(v, {valueMin: 0, valueMax: 100});
    expect(v.props['valueMin'], 0);
    expect(v.props['valueMax'], 100);
  });
});
