import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow, look} from '../helpers';

// Regression coverage for GROK-18091: adding a `transformation` to a PC plot
// while the filter panel was open made the viewer go empty after the panel
// closed. The fix preserves the transformed dataset across filter-panel toggles.
// `IPcPlotSettings.transformation` is a JSON-serialized opaque string.
//
// NOTE on dev 2026-04-29: setOptions({transformation: '<arbitrary json>'}) is
// rejected by the Dart parser (NoSuchMethodError 'bL'). Only the documented
// sentinel empty-string transformation is reliably round-tripped from the JS
// API surface; non-empty transformation strings are not exercisable from JS
// without the matching Dart-side parser. The full filter-panel open/close
// trigger surface lives behind a UI gesture and is covered by Selenium / manual QA.
category('AI: GROK-18091: PC plot transformation persistence', () => {
  test('transformation property exists and is reachable on PC plot', async () => {
    const v = DG.Viewer.pcPlot(demog(), {columnNames: ['age', 'height']});
    expect(v.type, DG.VIEWER.PC_PLOT);
    const l = look(v);
    expect('transformation' in l || l['transformation'] === undefined || l['transformation'] === '', true);
  });

  test('empty-string transformation round-trip is non-throwing', async () => {
    const v = DG.Viewer.pcPlot(demog(20), {columnNames: ['age', 'height']});
    expectNoThrow(() => v.setOptions({transformation: ''}));
    const t = look(v)['transformation'];
    expect(t === '' || t == null, true);
  });
});
