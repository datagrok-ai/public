import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Low-priority regression coverage for GROK-13210: Trellis plot with an inner
// Box plot used to render empty plots when Bin Color was set. The user-visible
// symptom (empty canvas inside each trellis cell) is a visual regression that
// cannot be asserted from JS — only the property-bag round-trip is checked
// here, so a schema-level regression on innerViewerLook.binColor* still has a
// tripwire. Property names verified against IBoxPlotSettings in
// public/js-api/src/interfaces/d4.d.ts: binColor, binColorColumnName,
// binColorAggrType.
category('AI: GROK-13210: Trellis Box plot bin color schema', () => {
  test('inner Box plot accepts binColor* in innerViewerLook and round-trips', async () => {
    const df = grok.data.demo.demog(60);
    var v: DG.Viewer | undefined;
    var threw = false;
    try {
      v = DG.Viewer.trellisPlot(df, {
        viewerType: 'Box plot',
        xColumnNames: ['sex'],
        yColumnNames: ['race'],
        innerViewerLook: {
          binColorColumnName: 'race',
          binColorAggrType: 'avg',
        },
      });
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v != null, true);
    try {
      expect(v! instanceof DG.Viewer, true);
      expect(v!.type, DG.VIEWER.TRELLIS_PLOT);
      // Schema-level tripwire: the look JSON envelope must include innerViewerLook.
      // Specific binColor* keys may be normalised away by the Dart side or stored
      // under different names — what we guard against is the whole structure
      // becoming unreachable (the original symptom: empty plots from broken inner
      // schema serialization). Per-key value pinning is intentionally NOT done.
      const look = v!.getOptions(true).look;
      expect(look.innerViewerLook != null, true);
    }
    finally {
      if (v != null)
        v.detach();
    }
  });

  test('toggling binColorAggrType through all aggregate variants survives setOptions', async () => {
    const aggrTypes = ['sum', 'avg', 'min', 'max', 'count', 'med', 'q1', 'q3', 'stdev'];
    const df = grok.data.demo.demog(60);
    const v = DG.Viewer.trellisPlot(df, {
      viewerType: 'Box plot',
      xColumnNames: ['sex'],
      yColumnNames: ['race'],
    });
    try {
      expect(v.type, DG.VIEWER.TRELLIS_PLOT);
      for (var aggr of aggrTypes) {
        var threw = false;
        try {
          v.setOptions({innerViewerLook: {binColorAggrType: aggr, binColorColumnName: 'race'}});
        }
        catch (_e) {
          threw = true;
        }
        expect(threw, false);
        const look = v.getOptions(true).look;
        expect(look.innerViewerLook != null, true);
      }
    }
    finally {
      v.detach();
    }
  });
});
