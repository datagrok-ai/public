import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectLook, expectNoThrow, look} from '../helpers';

// Regression coverage for GROK-18867: PC plot colorMin/colorMax persistence when columns change.
category('AI: GROK-18867: PC plot colorMin/colorMax persistence', () => {
  function coherent(a: any, b: any, av: number, bv: number): boolean {
    const bothNullish = a == null && b == null;
    const bothEqual = typeof a === 'number' && typeof b === 'number' && a === av && b === bv;
    return bothNullish || bothEqual;
  }

  const v = (): DG.PcPlot => DG.Viewer.pcPlot(demog(), {columnNames: ['age', 'height', 'weight']}) as DG.PcPlot;

  test('colorMin/colorMax round-trip into look after setOptions', async () => {
    const c = v();
    c.setOptions({colorColumnName: 'age', colorMin: 10, colorMax: 90});
    expectLook(c, {colorColumnName: 'age', colorMin: 10, colorMax: 90});
  });

  test('columnNames change preserves both color bounds coherently', async () => {
    const c = v();
    c.setOptions({colorColumnName: 'age', colorMin: 10, colorMax: 90});
    expectLook(c, {colorMin: 10, colorMax: 90});
    c.setOptions({columnNames: ['age', 'height']});
    const after = look(c);
    expect(coherent(after['colorMin'], after['colorMax'], 10, 90), true);
  });

  test('changing columnNames while color bounds are set does not throw', async () => {
    const c = v();
    c.setOptions({colorColumnName: 'age', colorMin: 10, colorMax: 90});
    expectNoThrow(() => c.setOptions({columnNames: ['height']}));
    const after = look(c);
    expect(coherent(after['colorMin'], after['colorMax'], 10, 90), true);
  });
});
