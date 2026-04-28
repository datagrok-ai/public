import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-17118 (1.23.0): density plot binShape +
// axis column type compatibility. With binShape='rectangle', switching
// xColumnName to a string column previously threw. The fix makes the
// viewer either ignore the invalid assignment or coerce it without
// crashing. These tests pin:
//   - binShape='rectangle' is accepted via setOptions and round-trips;
//   - changing xColumnName to a string column from rectangle state does
//     not throw, the viewer survives, and look.xColumnName ends up either
//     unchanged (rejected) or set to the string column (accepted);
//   - binShape round-trips for at least one non-rectangle valid value,
//     discovered from getProperties().choices when available, with a
//     fallback to known values like 'hex';
//   - getProperties() lists binShape, xColumnName, yColumnName on a
//     best-effort basis.
//
// Notes:
//   - viewer.dataFrame === df is unreliable, so we read rowCount from the
//     viewer for the survival check rather than identity-comparing dfs.
//   - binShape valid values aren't in a TS enum — we discover them at
//     runtime via property choices, with a static fallback.
category('AI: GROK-17118: Density plot binShape + axis column type', () => {
  test('binShape="rectangle" round-trips through setOptions/look', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.DENSITY_PLOT);

    v.setOptions({binShape: 'rectangle'});
    const look = v.getOptions(true).look;
    expect(look['binShape'], 'rectangle');
  });

  test('switching xColumnName to a string column from rectangle does not throw', async () => {
    const df = grok.data.demo.demog(50);
    // demog has numeric age/height and string columns like race/sex/site.
    const stringCol = 'race';
    expect(df.col(stringCol) != null, true);
    expect(df.col(stringCol)!.type !== DG.COLUMN_TYPE.INT &&
      df.col(stringCol)!.type !== DG.COLUMN_TYPE.FLOAT, true);

    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    v.setOptions({binShape: 'rectangle'});
    const beforeX = v.getOptions(true).look['xColumnName'];
    const beforeRowCount = v.dataFrame != null ? v.dataFrame.rowCount : df.rowCount;

    var threw = false;
    try {
      v.setOptions({xColumnName: stringCol});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    // Viewer instance must survive the assignment.
    expect(v instanceof DG.Viewer, true);
    const afterRowCount = v.dataFrame != null ? v.dataFrame.rowCount : df.rowCount;
    expect(afterRowCount, beforeRowCount);

    // The fix allows two valid outcomes: either retain the previous
    // numerical x column (rejected) or accept the string column. Either
    // way we must not crash and look must remain readable.
    const afterX = v.getOptions(true).look['xColumnName'];
    expect(afterX === beforeX || afterX === stringCol, true);
  });

  test('binShape round-trips for a non-rectangle valid value', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});

    // Discover valid binShape values via getProperties().choices when
    // available; otherwise fall back to known values.
    const props = v.getProperties() as DG.Property[];
    var choices: string[] | null = null;
    for (var p of props) {
      if (p.name === 'binShape') {
        try {
          const c = p.choices;
          if (Array.isArray(c) && c.length > 0)
            choices = c;
        }
        catch (_e) { /* property may not expose choices */ }
        break;
      }
    }

    const fallback = ['hex', 'rectangle'];
    const candidates = choices != null ? choices : fallback;
    var picked: string | null = null;
    for (var c of candidates) {
      if (c !== 'rectangle') { picked = c; break; }
    }
    if (picked == null) picked = 'hex';

    v.setOptions({binShape: picked});
    expect(v.getOptions(true).look['binShape'], picked);

    // And back to rectangle, to make sure switching between the two
    // valid values works without throwing.
    var threw = false;
    try {
      v.setOptions({binShape: 'rectangle'});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    expect(v.getOptions(true).look['binShape'], 'rectangle');
  });

  test('getProperties lists binShape/xColumnName/yColumnName (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.densityPlot(df, {x: 'age', y: 'height'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);

    const wanted = ['binShape', 'xColumnName', 'yColumnName'];
    const found: {[k: string]: boolean} = {};
    for (var p of props) {
      if (wanted.indexOf(p.name) >= 0)
        found[p.name] = true;
    }

    // Older descriptors may not list every prop; fall back to a look-bag
    // round-trip so the regression check is still meaningful.
    if (Object.keys(found).length < wanted.length) {
      v.setOptions({binShape: 'rectangle', xColumnName: 'age', yColumnName: 'height'});
      const look = v.getOptions(true).look;
      expect(look['binShape'], 'rectangle');
      expect(look['xColumnName'], 'age');
      expect(look['yColumnName'], 'height');
    }
    else {
      expect(found['binShape'], true);
      expect(found['xColumnName'], true);
      expect(found['yColumnName'], true);
    }
  });
});
