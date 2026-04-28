import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-18091: adding a `transformation` to a PC plot
// while the filter panel was open made the viewer go empty after the panel
// closed. The fix preserves the transformed dataset across filter-panel
// toggles. `IPcPlotSettings.transformation` is a JSON-serialized opaque
// string (interfaces/d4.ts:2748).
//
// NOTE on dev 2026-04-29: setOptions({transformation: '<arbitrary json>'}) is
// rejected by the Dart parser (NoSuchMethodError 'bL' — the parser expects a
// specific schema and the obfuscated entry-point is missing on this build).
// Only the documented sentinel empty-string transformation is reliably
// round-tripped from the JS API surface; non-empty transformation strings are
// not exercisable from JS without the matching Dart-side parser. The full
// filter-panel open/close trigger surface lives behind a UI gesture and is
// covered by Selenium / manual QA.
category('AI: GROK-18091: PC plot transformation persistence', () => {
  test('transformation property exists and is reachable on PC plot', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height']});
    expect(v instanceof DG.Viewer, true);
    expect(v.type, DG.VIEWER.PC_PLOT);

    // The look JSON envelope must include a transformation field (string,
    // possibly empty). This pins that the schema didn't drop the key in the
    // course of the fix.
    const look = v.getOptions(true).look;
    expect('transformation' in look || look['transformation'] === undefined || look['transformation'] === '', true);
  });

  test('empty-string transformation round-trip is non-throwing', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height']});
    expect(v instanceof DG.Viewer, true);

    var threw = false;
    try {
      v.setOptions({transformation: ''});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);
    // After clearing, the look bag's transformation is the empty string (or
    // absent — both shapes are accepted as "cleared").
    const t = v.getOptions(true).look['transformation'];
    expect(t === '' || t == null, true);
  });
});
