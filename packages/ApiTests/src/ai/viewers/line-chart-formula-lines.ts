import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectCleared, expectRoundTrip, look} from '../helpers';

// DG.LineChartViewer formula-lines reach — covered entirely through the
// base Viewer.meta.formulaLines helper (viewer.ts: ViewerFormulaLinesHelper)
// and the `formulaLines` / `showDataframeFormulaLines` / `showViewerFormulaLines`
// @Props on core/client/d4/lib/src/viewers/line_chart/line_chart_look.dart.
// No LineChart-specific Dart wraps: meta.formulaLines is inherited from the base
// Viewer and reads/writes through props['formulaLines'], so the previous
// FormulaLines_* wraps were redundant and have been removed.
category('AI: Viewers: LineChart Formula Lines', () => {
  const lineFormula = {title: 'AI test line', formula: '${height} = ${age} + 100', type: 'line'};
  const bandFormula = {title: 'AI test band', formula: '${height} in (160, 180)', type: 'band'};

  test('df.plot.line shorthand exposes meta.formulaLines helper', async () => {
    const df = demog();
    const v = df.plot.line({x: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
    expect(v instanceof DG.LineChartViewer, true);
    const fl = v.meta.formulaLines;
    expect(fl.items.length, 0);
    fl.add(lineFormula);
    expect(fl.items.length, 1);
    expect(fl.items[0].title, 'AI test line');
    // The helper persists straight into props['formulaLines'] (a @Prop).
    const parsed = JSON.parse(v.props['formulaLines']);
    expect(parsed.length, 1);
    expect(parsed[0].title, 'AI test line');
    fl.clear();
  });

  test('meta.formulaLines round-trip (add/items/clear) on LineChart', async () => {
    const v = DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height']});
    const fl = v.meta.formulaLines;
    fl.add(lineFormula);
    expect(fl.items.length, 1);
    fl.clear();
    expect(fl.items.length, 0);
  });

  test('meta.formulaLines.addBand keeps a band entry', async () => {
    const v = DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height']});
    const fl = v.meta.formulaLines;
    fl.addBand(bandFormula);
    expect(fl.items.length, 1);
    expect(fl.items[0].type, 'band');
    expect(fl.items[0].title, 'AI test band');
    fl.clear();
  });

  test('formulaLines round-trip via setOptions JSON survives + reflects in helper', async () => {
    const v = DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height']});
    const json = JSON.stringify([lineFormula]);
    expectRoundTrip(v, {formulaLines: json});
    expect(v.meta.formulaLines.items.length, 1);
    expect(v.meta.formulaLines.items[0].title, 'AI test line');
    // props mirror the round-tripped option.
    expect(v.props['formulaLines'], json);
    v.setOptions({formulaLines: ''});
    expectCleared(look(v)['formulaLines']);
  });

  test('showDataframeFormulaLines / showViewerFormulaLines are props that round-trip', async () => {
    const v = DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height']});
    expectRoundTrip(v, {showViewerFormulaLines: false});
    expectRoundTrip(v, {showViewerFormulaLines: true});
    expectRoundTrip(v, {showDataframeFormulaLines: false});
    expectRoundTrip(v, {showDataframeFormulaLines: true});
  });

  test('boundary: clear is idempotent + props stay empty', async () => {
    const v = DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height']});
    const fl = v.meta.formulaLines;
    fl.clear();
    expect(fl.items.length, 0);
    const raw = v.props['formulaLines'];
    expect(raw == null || raw === '' || raw === '[]', true);
    fl.clear();
    expect(fl.items.length, 0);
  });
}, {owner: 'agolovko@datagrok.ai'});
