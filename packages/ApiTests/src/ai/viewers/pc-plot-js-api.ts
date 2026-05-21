import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectPropAndLook, expectRoundTrip, findProp, look, subscribeAll, withAttachedViewer} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:275 (DG.Viewer.pcPlot),
// public/js-api/src/viewer.ts:704 (DG.PcPlot — bare class, no Viewer suffix),
// public/js-api/src/interfaces/d4.ts:2692 (IPcPlotSettings).
// PC-only JS surface: factory DG.Viewer.pcPlot, the defining
// `columnNames: Array<string>` axis-list round-trip via setOptions and look,
// PC-specific Boolean toggles (normalizeEachColumn / showCurrentLine /
// showAllLines), the PC-only event Observables (onLineClicked / onLineHovered),
// and view.addViewer(VIEWER.PC_PLOT) returning a typed DG.PcPlot. Note: there
// is no df.plot.pcPlot shorthand on DataFramePlotHelper.
category('AI: Viewers: PC Plot JS API', () => {
  test('factory DG.Viewer.pcPlot returns typed DG.PcPlot with PC settings', async () => {
    const df = demog();
    const v = DG.Viewer.pcPlot(df, {columnNames: ['age', 'height', 'weight']});
    expect(v instanceof DG.PcPlot, true);
    expect(v.dataFrame === df, true);
    expectArray(v.props['columnNames'], ['age', 'height', 'weight']);
  });

  test('columnNames array round-trip via getOptions and setOptions', async () => {
    const v = DG.Viewer.pcPlot(demog(), {columnNames: ['age', 'height']});
    expectArray(v.props['columnNames'], ['age', 'height']);
    expectArray(look(v)['columnNames'] as string[], ['age', 'height']);
    expectRoundTrip(v, {columnNames: ['age', 'weight']});
  });

  test('PC-specific Boolean props round-trip and getProperties surfaces columnNames', async () => {
    const v = DG.Viewer.pcPlot(demog(), {columnNames: ['age', 'height', 'weight']});
    const flipped = {
      normalizeEachColumn: !v.props['normalizeEachColumn'],
      showCurrentLine: !v.props['showCurrentLine'],
      showAllLines: !v.props['showAllLines'],
    };
    v.setOptions(flipped);
    expectPropAndLook(v, flipped);
    expect(findProp(v, 'columnNames') != null, true);
  });

  test('onLineClicked and onLineHovered are rxjs Observables', async () => {
    const v = DG.Viewer.pcPlot(demog(20), {columnNames: ['age', 'height']}) as DG.PcPlot;
    subscribeAll([v.onLineClicked, v.onLineHovered])();
  });

  test('view.addViewer attaches a typed DG.PcPlot', async () => {
    await withAttachedViewer<DG.PcPlot>(demog(), DG.VIEWER.PC_PLOT, {columnNames: ['age', 'height', 'weight']},
      (v) => expect(v instanceof DG.PcPlot, true));
  });
});
