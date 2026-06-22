import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectCleared, expectRoundTrip, look, until, withTableView} from '../helpers';

// LineChart prop JSON-shape round-trips not covered by line-chart-js-api.ts.
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

  // whiskersType flip resets markerOpacity to its optimal value.
  test('whiskersType change auto-adjusts markerOpacity', async () => {
    await withTableView(demog(), async (tv) => {
      const c = v();
      tv.addViewer(c);
      c.setOptions({whiskersType: 'None'});
      await until(() => look(c)['markerOpacity'] === 100);
      expect(look(c)['markerOpacity'], 100);
      c.setOptions({whiskersType: 'Med | Q1, Q3'});
      await until(() => look(c)['markerOpacity'] !== 100);
      expect(look(c)['markerOpacity'] !== 100, true);
      c.setOptions({whiskersType: 'None'});
      await until(() => look(c)['markerOpacity'] === 100);
      expect(look(c)['markerOpacity'], 100);
    });
  });

  // Optimal markerSize needs a real layout pass, so assert direction (user value overridden), not the exact number.
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
