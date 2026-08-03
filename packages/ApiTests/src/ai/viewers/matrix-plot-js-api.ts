import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {demog, df, look, expectRoundTripPropAndLook, expectAttachesNoThrow} from '../helpers';

// DG.MatrixPlot

category('AI: Viewers: Matrix Plot', () => {
  test('matrixPlot factory returns typed DG.MatrixPlot of the right type', async () => {
    const v = DG.Viewer.matrixPlot(demog(), {xColumnNames: ['age'], yColumnNames: ['height']});
    expect(v instanceof DG.MatrixPlot, true);
    expect(v.type, DG.VIEWER.MATRIX_PLOT);
  });

  test('cellPlotType + style props round-trip', async () => {
    const v = DG.Viewer.matrixPlot(demog());
    // cellPlotType valid values are 'Density plot' (default) and 'Scatter plot' (@Prop choices).
    expectRoundTripPropAndLook(v, {cellPlotType: 'Scatter plot', autoLayout: false, showXAxes: true, showYAxes: true});
  });

  test('x/y column array props round-trip via setOptions + expectArray', async () => {
    const v = DG.Viewer.matrixPlot(demog());
    // xColumnNames/yColumnNames are List<String> — compare with a deep array check,
    // not the scalar expect() inside expectPropAndLook.
    v.setOptions({xColumnNames: ['age', 'height'], yColumnNames: ['weight']});
    expectArray(look(v).xColumnNames, ['age', 'height']);
    expectArray(look(v).yColumnNames, ['weight']);
  });

  // onViewerRendered is a base-Viewer event (covered in Viewers: Lifecycle Events); the boundary
  // test below covers attach. No per-viewer-type render-event re-test here.

  test('boundary: single-row frame and a minimal no-explicit-columns frame both attach without throwing', async () => {
    const one = df([['a', DG.COLUMN_TYPE.FLOAT, [1.0]], ['b', DG.COLUMN_TYPE.FLOAT, [2.0]]]);
    await expectAttachesNoThrow(one, DG.VIEWER.MATRIX_PLOT, {xColumnNames: ['a'], yColumnNames: ['b']});
    const numeric = df([['x', DG.COLUMN_TYPE.FLOAT, [1, 2, 3]], ['y', DG.COLUMN_TYPE.FLOAT, [4, 5, 6]]]);
    await expectAttachesNoThrow(numeric, DG.VIEWER.MATRIX_PLOT, {});
  });
}, {owner: 'agolovko@datagrok.ai'});
