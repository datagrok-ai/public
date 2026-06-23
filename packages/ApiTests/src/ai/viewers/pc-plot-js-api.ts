import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, expectPropAndLook, expectRoundTrip, findProp, look, subscribeAll, withAttachedViewer} from '../helpers';

// PcPlot JS API: factory, columnNames round-trip, bool toggles, events, addViewer.
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
