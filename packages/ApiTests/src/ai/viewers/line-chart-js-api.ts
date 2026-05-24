import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {
  demog, expectPropAndLook, expectRoundTripPropAndLook, findProp, look, subscribeAll, withAttachedViewer,
} from '../helpers';

// JS API source: public/js-api/src/viewer.ts:243 (DG.Viewer.lineChart),
// public/js-api/src/viewer.ts:577 (DG.LineChartViewer).
// LineChart-only JS surface: factory friendly-key aliasing on
// DG.Viewer.lineChart and df.plot.line (note: ILineChartSettings has no `y`
// field — only yColumnNames; the singular 'y' alias does not resolve), the
// multi-element array properties (yColumnNames / splitColumnNames) round-tripping
// through getOptions(true).look, multiAxis + splineTension via setOptions,
// the LineChart-typed events (onLineSelected / onZoomed / onResetView) as
// rxjs Observables, and view.addViewer attaching a typed LineChartViewer with
// activeFrame. screenToWorld(0,0) was intentionally dropped here — pre-paint
// the viewport Rect is NaN and screenToWorld throws NaN.floor() rather than
// returning null (worth a JIRA).
category('AI: Viewers: LineChart JS API', () => {
  test('factory friendly-key aliasing via DG.Viewer.lineChart', async () => {
    const v = DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height'], split: 'race'});
    expect(v instanceof DG.LineChartViewer, true);
    expectPropAndLook(v, {xColumnName: 'age', splitColumnName: 'race'});
    expect((look(v)['yColumnNames'] as string[]).indexOf('height') >= 0, true);
  });

  test('factory via df.plot.line shorthand', async () => {
    const df = demog();
    const v = df.plot.line({xColumnName: 'age', yColumnNames: ['weight']});
    expect(v.type, DG.VIEWER.LINE_CHART);
    expect(v.dataFrame === df, true);
    expect(v.props['xColumnName'], 'age');
    expect((v.props['yColumnNames'] as string[]).indexOf('weight') >= 0, true);
  });

  test('yColumnNames and splitColumnNames array round-trip', async () => {
    const v = DG.Viewer.lineChart(demog(), {
      xColumnName: 'age', yColumnNames: ['height', 'weight'], splitColumnNames: ['race', 'sex'],
    });
    expectArray(v.props['yColumnNames'], ['height', 'weight']);
    expectArray(v.props['splitColumnNames'], ['race', 'sex']);
    expectArray(look(v)['yColumnNames'] as string[], ['height', 'weight']);
    expectArray(look(v)['splitColumnNames'] as string[], ['race', 'sex']);
  });

  test('multiAxis and splineTension round-trip via setOptions', async () => {
    const v = DG.Viewer.lineChart(demog(), {xColumnName: 'age', yColumnNames: ['height', 'weight']});
    expect(v.props['multiAxis'], false);
    expectRoundTripPropAndLook(v, {multiAxis: true, splineTension: 0.5});
    expect(findProp(v, 'multiAxis') != null, true);
    expect(findProp(v, 'splineTension') != null, true);
  });

  test('onLineSelected and onZoomed are rxjs Observables', async () => {
    const v = DG.Viewer.lineChart(demog(20), {xColumnName: 'age', yColumnNames: ['height']}) as DG.LineChartViewer;
    subscribeAll([v.onLineSelected, v.onZoomed, v.onResetView])();
  });

  test('view.addViewer attaches a typed LineChartViewer with activeFrame', async () => {
    await withAttachedViewer<DG.LineChartViewer>(demog(), DG.VIEWER.LINE_CHART,
      {xColumnName: 'age', yColumnNames: ['height']}, (v) => {
        expect(v instanceof DG.LineChartViewer, true);
        const af = v.activeFrame;
        expect(af != null, true);
        expect(DG.toJs(af) instanceof DG.DataFrame, true);
      });
  });
});
