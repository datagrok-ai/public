import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-19760 (also covers GROK-19581 — same root cause).
// Setting an inverted, NaN, or zero-range valueMin/valueMax used to crash the
// histogram. The fix is defensive: the viewer must survive, with no exception
// escaping setOptions. The exact clamp / rejection policy is intentionally not
// pinned — we only assert no-throw, viewer liveness (via rowCount), and
// props↔look round-trip.
category('AI: GROK-19760: Histogram invalid valueMin/valueMax', () => {
  test('inverted range does not throw and survives', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.histogram({value: 'age'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.HISTOGRAM);

    var threw = false;
    try {
      v.setOptions({valueMin: 100, valueMax: 0});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    // Liveness: viewer is still attached to a df with the same rowCount.
    // (We compare rowCount rather than identity — the dataFrame getter can
    // return a fresh wrapper across the Dart→JS boundary.)
    expect(v.dataFrame != null, true);
    expect(v.dataFrame!.rowCount, df.rowCount);
    // Round-trip: look bag agrees with props for both keys.
    const look = v.getOptions(true).look;
    expect(v.props['valueMin'], look['valueMin']);
    expect(v.props['valueMax'], look['valueMax']);
  });

  test('NaN valueMin does not throw', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.histogram({value: 'age'});
    var threw = false;
    try {
      v.setOptions({valueMin: NaN, valueMax: 50});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.dataFrame != null, true);
    expect(v.dataFrame!.rowCount, df.rowCount);
    // Round-trip the look bag against props — values may be clamped or
    // rejected (null/undefined). Normalise both sides since props returns null
    // and look returns undefined for cleared values.
    const look = v.getOptions(true).look;
    expect(v.props['valueMin'] == look['valueMin'] || (v.props['valueMin'] == null && look['valueMin'] == null), true);
    expect(v.props['valueMax'] == look['valueMax'] || (v.props['valueMax'] == null && look['valueMax'] == null), true);
  });

  test('zero range survives, then a valid range lands cleanly', async () => {
    const df = grok.data.demo.demog(50);
    const v = df.plot.histogram({value: 'age'});

    var threw = false;
    try {
      v.setOptions({valueMin: 50, valueMax: 50});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.dataFrame != null, true);
    expect(v.dataFrame!.rowCount, df.rowCount);

    // Now a valid range — must land cleanly in look and round-trip with props.
    threw = false;
    try {
      v.setOptions({valueMin: 0, valueMax: 100});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    const look = v.getOptions(true).look;
    expect(look['valueMin'], 0);
    expect(look['valueMax'], 100);
    expect(v.props['valueMin'], 0);
    expect(v.props['valueMax'], 100);
  });
});
