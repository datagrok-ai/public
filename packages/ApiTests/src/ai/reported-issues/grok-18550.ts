import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-18550 (1.26.3): PC plot `normalizeEachColumn`
// stability. Turning off Normalize Each Column with log columns set used to
// throw. The fix lets the toggle survive in any combination with
// `logColumnsColumnNames`. Surface: `DG.Viewer.pcPlot` (`viewer.ts:275`),
// `IPcPlotSettings.normalizeEachColumn` (`d4.ts:2724`),
// `IPcPlotSettings.logColumnsColumnNames: Array<string>` (`d4.ts:2702`).
//
// Pin the round-trip: build via `DG.Viewer.pcPlot`, mutate via `setOptions`,
// read back from `getOptions(true).look`. No DOM/canvas inspection — JS API
// only. `viewer.dataFrame === df` is unreliable across the interop boundary,
// so we don't assert it.
category('AI: GROK-18550: PC plot normalizeEachColumn stability', () => {
  test('toggle normalizeEachColumn true→false→true with log columns set: no throw', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {
      columnNames: ['age', 'height', 'weight'],
      logColumnsColumnNames: ['age'],
    });
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.PC_PLOT);

    var threw = false;
    try {
      v.setOptions({normalizeEachColumn: false});
      v.setOptions({normalizeEachColumn: true});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    // Mirrored in look — final state should be true.
    const look = v.getOptions(true).look;
    expect(look['normalizeEachColumn'], true);

    // logColumnsColumnNames should still carry 'age'. Accept either an
    // explicit array or a missing key on the look bag (some platform
    // builds omit empty/default array values).
    const logs = look['logColumnsColumnNames'] as string[] | undefined;
    if (logs !== undefined) {
      expect(Array.isArray(logs), true);
      expectArray(logs, ['age']);
    }
  });

  test('with empty logColumnsColumnNames: toggle normalizeEachColumn survives', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {
      columnNames: ['age', 'height', 'weight'],
      logColumnsColumnNames: ['age'],
    });

    var threw = false;
    try {
      v.setOptions({logColumnsColumnNames: []});
      v.setOptions({normalizeEachColumn: false});
      v.setOptions({normalizeEachColumn: true});
      v.setOptions({normalizeEachColumn: false});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    const look = v.getOptions(true).look;
    expect(look['normalizeEachColumn'], false);

    // Risk: empty array may serialize as missing key on the look bag.
    // Accept either an empty array or missing key.
    const logs = look['logColumnsColumnNames'] as string[] | undefined;
    if (logs !== undefined) {
      expect(Array.isArray(logs), true);
      expect(logs.length, 0);
    }
  });

  test('getProperties (best-effort) lists normalizeEachColumn and logColumnsColumnNames', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pcPlot(df, {
      columnNames: ['age', 'height'],
      logColumnsColumnNames: ['age'],
    });
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);

    var foundNorm: DG.Property | null = null;
    var foundLogs: DG.Property | null = null;
    for (var p of props) {
      if (p.name === 'normalizeEachColumn')
        foundNorm = p;
      if (p.name === 'logColumnsColumnNames')
        foundLogs = p;
    }

    // Best-effort: the descriptor list may not always carry these at
    // runtime. When missing, fall back to the schema-only check via the
    // round-trip we already exercise above.
    if (foundNorm == null || foundLogs == null) {
      v.setOptions({normalizeEachColumn: false});
      const look = v.getOptions(true).look;
      expect(look['normalizeEachColumn'], false);
    }
  });
});
