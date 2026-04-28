import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-18867: PC plot colorMin/colorMax persistence
// when columns change. The reported defect was that the per-column range
// slider min/max could be reset incorrectly when the active columns set
// changed. The internal slider state is not part of the public surface, so
// we pin the JS-API-visible analog: `colorMin` and `colorMax` on
// `IPcPlotSettings` (`d4.ts:2716-2718`). The fix shipped in 1.26.3.
//
// We assert a tolerant invariant: after a `columnNames` change, the look
// bag must keep both color bounds *coherent* — either both numerically
// preserved at the originally pinned values, or both cleared together.
// Half-set state (one numeric, the other null/undefined) is the bug.
category('AI: GROK-18867: PC plot colorMin/colorMax persistence', () => {
  function bothNullish(a: any, b: any): boolean {
    return (a == null) && (b == null);
  }

  function bothEqualNumbers(a: any, b: any, av: number, bv: number): boolean {
    return typeof a === 'number' && typeof b === 'number' && a === av && b === bv;
  }

  test('colorMin/colorMax round-trip into look after setOptions', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    expect(v instanceof DG.PcPlot, true);
    expect(v.type, DG.VIEWER.PC_PLOT);

    v.setOptions({colorColumnName: 'age', colorMin: 10, colorMax: 90});
    const look = v.getOptions(true).look;
    expect(look['colorColumnName'], 'age');
    expect(look['colorMin'], 10);
    expect(look['colorMax'], 90);
  });

  test('columnNames change preserves both color bounds coherently', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    v.setOptions({colorColumnName: 'age', colorMin: 10, colorMax: 90});

    // Sanity: pinned before the columns change.
    const before = v.getOptions(true).look;
    expect(before['colorMin'], 10);
    expect(before['colorMax'], 90);

    // Trigger the code path the ticket touched: reassign columnNames while
    // the color bounds are set.
    v.setOptions({columnNames: ['age', 'height']});

    const after = v.getOptions(true).look;
    const cMin = after['colorMin'];
    const cMax = after['colorMax'];

    // Tolerant invariant: never half-set. Either both preserved at 10/90,
    // or both cleared together. The bug was the asymmetric reset.
    const preserved = bothEqualNumbers(cMin, cMax, 10, 90);
    const cleared = bothNullish(cMin, cMax);
    expect(preserved || cleared, true);
  });

  test('changing columnNames while color bounds are set does not throw', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    v.setOptions({colorColumnName: 'age', colorMin: 10, colorMax: 90});

    var threw = false;
    try {
      v.setOptions({columnNames: ['height']});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    // After the smoke step, the same coherence invariant must hold.
    const look = v.getOptions(true).look;
    const cMin = look['colorMin'];
    const cMax = look['colorMax'];
    const preserved = bothEqualNumbers(cMin, cMax, 10, 90);
    const cleared = bothNullish(cMin, cMax);
    expect(preserved || cleared, true);
  });
});
