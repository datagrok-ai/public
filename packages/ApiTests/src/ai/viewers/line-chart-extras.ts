import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectCleared, expectRoundTrip, look, until, wait, withTableView} from '../helpers';

// JS API source: public/js-api/src/viewer.ts (DG.Viewer.lineChart),
// public/js-api/src/interfaces/d4.d.ts:686 (ILineChartSettings),
// core/client/d4/lib/src/viewers/line_chart/line_chart_look.dart (LineChartLook).
// Line-chart-only JSON-shape coverage for props the existing
// src/ai/viewers/line-chart-js-api.ts (factory aliases, multi-Y arrays,
// multiAxis + splineTension, events, addViewer + activeFrame) does not pin:
// xAxisType/yAxisType (AxisType linear/logarithmic), xAxisLabelOrientation
// choices ('Auto'/'Horz'/'Vert'), the invert/show* axis-and-grid bool
// envelope, yAxisTitle/y2AxisTitle string round-trip, markerOpacity 0..100
// boundary, packCategories/yGlobalScale/axesFollowFilter/showCurrentRowLine/
// showMouseOverCategory combined bool envelope, and overviewColumnName +
// overviewAggrType combined round-trip. All assertions read state via
// getOptions(true).look — no first-paint geometry.
category('AI: Viewers: LineChart extras', () => {
  const v = (): DG.LineChartViewer => DG.Viewer.lineChart(demog(), {x: 'age', yColumnNames: ['height']});

  test('xAxisType + yAxisType (AxisType) round-trip + getProperties choices', async () => {
    const c = v();
    expectRoundTrip(c, {xAxisType: 'logarithmic', yAxisType: 'logarithmic'});
    expectRoundTrip(c, {xAxisType: 'linear', yAxisType: 'linear'});
    expectChoices(c, 'yAxisType', ['linear', 'logarithmic']);
  });

  test('xAxisLabelOrientation choices Auto/Horz/Vert round-trip + getProperties choices', async () => {
    const c = v();
    for (const o of ['Horz', 'Vert', 'Auto'])
      expectRoundTrip(c, {xAxisLabelOrientation: o});
    expectChoices(c, 'xAxisLabelOrientation', ['Auto', 'Horz', 'Vert']);
  });

  test('invertXAxis/show*Axis/show*GridLines combined bool round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {invertXAxis: true, showXAxis: false, showYAxis: false,
      showVerticalGridLines: false, showHorizontalGridLines: false});
    expectRoundTrip(c, {invertXAxis: false, showXAxis: true, showYAxis: true,
      showVerticalGridLines: true, showHorizontalGridLines: true});
  });

  test('yAxisTitle + y2AxisTitle string round-trip + clearing via empty string', async () => {
    const c = v();
    expectRoundTrip(c, {yAxisTitle: 'Height (cm)', y2AxisTitle: 'Weight (kg)'});
    c.setOptions({yAxisTitle: '', y2AxisTitle: ''});
    const cleared = look(c);
    expectCleared(cleared['yAxisTitle']);
    expectCleared(cleared['y2AxisTitle']);
  });

  test('markerOpacity 0..100 boundary round-trip', async () => {
    const c = v();
    for (const n of [0, 50, 100])
      expectRoundTrip(c, {markerOpacity: n});
  });

  test('packCategories/yGlobalScale/axesFollowFilter/showCurrentRowLine/showMouseOverCategory bools', async () => {
    const c = v();
    expectRoundTrip(c, {packCategories: false, yGlobalScale: true,
      axesFollowFilter: false, showCurrentRowLine: true, showMouseOverCategory: false});
    expectRoundTrip(c, {packCategories: true, yGlobalScale: false,
      axesFollowFilter: true, showCurrentRowLine: false, showMouseOverCategory: true});
  });

  test('overviewColumnName + overviewAggrType combined round-trip', async () => {
    const c = v();
    expectRoundTrip(c, {overviewColumnName: 'weight', overviewAggrType: 'med'});
    c.setOptions({overviewColumnName: '', overviewAggrType: 'avg'});
    const cleared = look(c);
    expectCleared(cleared['overviewColumnName']);
    expect(cleared['overviewAggrType'], 'avg');
  });

  // Property-dependency coverage (commit cd5cfcb994, "Line chart: Added
  // marker opacity adjusting based on whisker visibility"). When
  // `whiskersType` flips, line_chart_core.dart sets `markerOpacity` to its
  // optimal value (100 with no whiskers, 35 with whiskers).
  test('whiskersType change auto-adjusts markerOpacity', async () => {
    await withTableView(demog(), async (tv) => {
      const c = v();
      tv.addViewer(c);
      c.setOptions({whiskersType: 'None'});
      await wait(100);
      expect(look(c)['markerOpacity'], 100);
      c.setOptions({whiskersType: 'Med | Q1, Q3'});
      await wait(100);
      expect(look(c)['markerOpacity'] !== 100, true);
      c.setOptions({whiskersType: 'None'});
      await wait(100);
      expect(look(c)['markerOpacity'], 100);
    });
  });

  // Property-dependency coverage. When `autoLayout` flips to true and the
  // current `markerSize` differs from the chart's optimal size,
  // line_chart_core.dart resets `markerSize` to optimal. Optimal depends on
  // `chartsBox.area / dataFrame.rowCount` which is only valid after a layout
  // pass, so we drive this through a real TableView and assert the *direction*
  // (user-set value gets overridden), not the exact optimal number.
  test('autoLayout=true overrides a user-set markerSize', async () => {
    await withTableView(demog(), async (tv) => {
      const c = v();
      tv.addViewer(c);
      await until(() => c.root.querySelector('canvas') != null);
      c.setOptions({autoLayout: false, markerSize: 99});
      await until(() => look(c)['markerSize'] === 99);
      expect(look(c)['markerSize'], 99);
      c.setOptions({autoLayout: true});
      await until(() => look(c)['markerSize'] !== 99);
      const after = look(c);
      expect(after['autoLayout'], true);
      expect(after['markerSize'] !== 99 && after['markerSize'] >= 1 && after['markerSize'] <= 5, true);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
