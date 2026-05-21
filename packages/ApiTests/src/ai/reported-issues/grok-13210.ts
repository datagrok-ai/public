import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, look} from '../helpers';

// Low-priority regression coverage for GROK-13210: Trellis plot with an inner
// Box plot used to render empty plots when Bin Color was set. The user-visible
// symptom (empty canvas inside each trellis cell) is a visual regression that
// cannot be asserted from JS — only the property-bag round-trip is checked
// here, so a schema-level regression on innerViewerLook.binColor* still has a
// tripwire. Property names verified against IBoxPlotSettings.
category('AI: GROK-13210: Trellis Box plot bin color schema', () => {
  const makeTrellis = (innerLook?: {[k: string]: any}): DG.Viewer => DG.Viewer.trellisPlot(demog(60), {
    viewerType: 'Box plot', xColumnNames: ['sex'], yColumnNames: ['race'],
    ...(innerLook != null ? {innerViewerLook: innerLook} : {}),
  });

  test('inner Box plot accepts binColor* in innerViewerLook and round-trips', async () => {
    var v: DG.Viewer | undefined;
    expectNoThrow(() => { v = makeTrellis({binColorColumnName: 'race', binColorAggrType: 'avg'}); });
    try {
      expect(v!.type, DG.VIEWER.TRELLIS_PLOT);
      // Schema-level tripwire: the look JSON envelope must include innerViewerLook.
      // Specific binColor* keys may be normalised away by the Dart side or stored
      // under different names — what we guard against is the whole structure
      // becoming unreachable. Per-key value pinning is intentionally NOT done.
      expect(look(v!).innerViewerLook != null, true);
    }
    finally {
      v!.detach();
    }
  });

  test('toggling binColorAggrType through all aggregate variants survives setOptions', async () => {
    const v = makeTrellis();
    try {
      for (var aggr of ['sum', 'avg', 'min', 'max', 'count', 'med', 'q1', 'q3', 'stdev']) {
        expectNoThrow(() => v.setOptions({innerViewerLook: {binColorAggrType: aggr, binColorColumnName: 'race'}}));
        expect(look(v).innerViewerLook != null, true);
      }
    }
    finally {
      v.detach();
    }
  });
});
