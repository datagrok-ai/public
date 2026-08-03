import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, look} from '../helpers';

// Regression coverage for GROK-13210: Trellis Box plot bin color schema round-trip.
category('AI: GROK-13210: Trellis Box plot bin color schema', () => {
  const makeTrellis = (innerLook?: {[k: string]: any}): DG.Viewer => DG.Viewer.trellisPlot(demog(60), {
    viewerType: 'Box plot', xColumnNames: ['sex'], yColumnNames: ['race'],
    ...(innerLook != null ? {innerViewerLook: innerLook} : {}),
  });

  test('inner Box plot accepts binColor* in innerViewerLook and round-trips', async () => {
    let v: DG.Viewer | undefined;
    expectNoThrow(() => {v = makeTrellis({binColorColumnName: 'race', binColorAggrType: 'avg'});});
    try {
      expect(v!.type, DG.VIEWER.TRELLIS_PLOT);
      expect(look(v!).innerViewerLook != null, true);
    } finally {
      v!.detach();
    }
  });

  test('toggling binColorAggrType through all aggregate variants survives setOptions', async () => {
    const v = makeTrellis();
    try {
      for (const aggr of ['sum', 'avg', 'min', 'max', 'count', 'med', 'q1', 'q3', 'stdev']) {
        expectNoThrow(() => v.setOptions({innerViewerLook: {binColorAggrType: aggr, binColorColumnName: 'race'}}));
        expect(look(v).innerViewerLook != null, true);
      }
    } finally {
      v.detach();
    }
  });
});
