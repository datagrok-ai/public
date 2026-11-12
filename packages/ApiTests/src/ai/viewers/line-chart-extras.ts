import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectChoices, expectCleared, expectFiresWithin, expectRoundTrip, look, reloadLayoutWithLook,
  until, withTableView} from '../helpers';

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

  test('opacity and transparencies are saved as the color alpha', async () => {
    const c = v();
    c.setOptions({regressionLineTransparency: 0.2, movingAverageLineTransparency: 0.6});
    expect(DG.Color.a(c.props.regressionLineColor), 204);
    expect(DG.Color.a(c.props.movingAverageLineColor), 102);

    await withTableView(demog(), (tv) => {
      tv.lineChart({x: 'age', yColumnNames: ['height']});
      const lc = reloadLayoutWithLook(tv, DG.VIEWER.LINE_CHART,
        {regressionLineTransparency: 0.5, regressionLineColor: 0xFFFF0000, lineTransparency: 0.75, markerOpacity: 50});
      expect((lc.props.regressionLineColor >>> 0).toString(16), '80ff0000');
      expect(DG.Color.a(lc.props.lineColor), 64);
      expect(DG.Color.a(lc.props.markerColor), 127);
      const whiskered = reloadLayoutWithLook(tv, DG.VIEWER.LINE_CHART,
        {whiskersType: 'Med | Q1, Q3', markerOpacity: null});
      expect(DG.Color.a(whiskered.props.markerColor), 89);
    });
  });

  test('setMarkerOpacity sets the alpha of markerColor', async () => {
    const c = v();
    expect(c.getProperties().some((p) => p.name === 'markerOpacity' || p.name === 'lineTransparency'), false);
    c.setMarkerOpacity(50);
    expect(DG.Color.a(c.props.markerColor), 127);
    c.setMarkerOpacity(0);
    expect(DG.Color.a(c.props.markerColor), 1);
    expect(DG.Color.a(c.props.lineColor), 255);
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

  // whiskersType flip resets the alpha of markerColor to its optimal value.
  test('whiskersType change auto-adjusts the marker color alpha', async () => {
    await withTableView(demog(), async (tv) => {
      const c = v();
      tv.addViewer(c);
      await until(() => c.root.querySelector('canvas') != null);
      expect(DG.Color.a(c.props.markerColor), 255); // no whiskers -> full opacity
      await expectFiresWithin(c.onAfterDrawScene, () => c.setOptions({whiskersType: 'Med | Q1, Q3'}), 2000);
      expect(DG.Color.a(c.props.markerColor), 89);
      await expectFiresWithin(c.onAfterDrawScene, () => c.setOptions({whiskersType: 'None'}), 2000);
      expect(DG.Color.a(c.props.markerColor), 255);
    });
  });

  // Optimal markerSize needs a real layout pass, so assert direction (user value overridden), not the exact number.
  test('autoLayout=true overrides a user-set markerSize', async () => {
    await withTableView(demog(), async (tv) => {
      const c = v();
      tv.addViewer(c);
      await until(() => c.root.querySelector('canvas') != null);
      c.setOptions({autoLayout: false, markerSize: 25});
      await until(() => look(c)['markerSize'] === 25);
      expect(look(c)['markerSize'], 25);
      c.setOptions({autoLayout: true});
      await until(() => look(c)['markerSize'] !== 25);
      const after = look(c);
      expect(after['autoLayout'], true);
      expect(after['markerSize'] !== 25 && after['markerSize'] >= 1 && after['markerSize'] <= 5, true);
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
