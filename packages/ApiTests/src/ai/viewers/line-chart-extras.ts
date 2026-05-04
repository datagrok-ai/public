import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

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
// getOptions(true).look — no first-paint geometry (round 8's
// screenToWorld NaN.floor() pitfall stays clear).
category('AI: Viewers: LineChart extras', () => {
  // no negative case: LineChart setOptions has no defined failure mode for
  // these look props — invalid values are coerced or ignored, not thrown.

  test('xAxisType + yAxisType (AxisType) round-trip + getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({xAxisType: 'logarithmic', yAxisType: 'logarithmic'});
    const log = v.getOptions(true).look;
    expect(log['xAxisType'], 'logarithmic');
    expect(log['yAxisType'], 'logarithmic');
    v.setOptions({xAxisType: 'linear', yAxisType: 'linear'});
    const lin = v.getOptions(true).look;
    expect(lin['xAxisType'], 'linear');
    expect(lin['yAxisType'], 'linear');
    const props = v.getProperties();
    var yProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'yAxisType') { yProp = p; break; }
    expect(yProp != null, true);
    expect(Array.isArray(yProp!.choices), true);
    expect(yProp!.choices.indexOf('linear') >= 0, true);
    expect(yProp!.choices.indexOf('logarithmic') >= 0, true);
  });

  test('xAxisLabelOrientation choices Auto/Horz/Vert round-trip + getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({xAxisLabelOrientation: 'Horz'});
    expect(v.getOptions(true).look['xAxisLabelOrientation'], 'Horz');
    v.setOptions({xAxisLabelOrientation: 'Vert'});
    expect(v.getOptions(true).look['xAxisLabelOrientation'], 'Vert');
    v.setOptions({xAxisLabelOrientation: 'Auto'});
    expect(v.getOptions(true).look['xAxisLabelOrientation'], 'Auto');
    const props = v.getProperties();
    var orientProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'xAxisLabelOrientation') { orientProp = p; break; }
    expect(orientProp != null, true);
    expect(Array.isArray(orientProp!.choices), true);
    expect(orientProp!.choices.indexOf('Auto') >= 0, true);
    expect(orientProp!.choices.indexOf('Horz') >= 0, true);
    expect(orientProp!.choices.indexOf('Vert') >= 0, true);
  });

  test('invertXAxis/show*Axis/show*GridLines combined bool round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({
      invertXAxis: true, showXAxis: false, showYAxis: false,
      showVerticalGridLines: false, showHorizontalGridLines: false,
    });
    const off = v.getOptions(true).look;
    expect(off['invertXAxis'], true);
    expect(off['showXAxis'], false);
    expect(off['showYAxis'], false);
    expect(off['showVerticalGridLines'], false);
    expect(off['showHorizontalGridLines'], false);
    v.setOptions({
      invertXAxis: false, showXAxis: true, showYAxis: true,
      showVerticalGridLines: true, showHorizontalGridLines: true,
    });
    const on = v.getOptions(true).look;
    expect(on['invertXAxis'], false);
    expect(on['showXAxis'], true);
    expect(on['showYAxis'], true);
    expect(on['showVerticalGridLines'], true);
    expect(on['showHorizontalGridLines'], true);
  });

  test('yAxisTitle + y2AxisTitle string round-trip + clearing via empty string', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({yAxisTitle: 'Height (cm)', y2AxisTitle: 'Weight (kg)'});
    const set = v.getOptions(true).look;
    expect(set['yAxisTitle'], 'Height (cm)');
    expect(set['y2AxisTitle'], 'Weight (kg)');
    v.setOptions({yAxisTitle: '', y2AxisTitle: ''});
    const cleared = v.getOptions(true).look;
    expect(cleared['yAxisTitle'] === '' || cleared['yAxisTitle'] == null, true);
    expect(cleared['y2AxisTitle'] === '' || cleared['y2AxisTitle'] == null, true);
  });

  test('markerOpacity 0..100 boundary round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({markerOpacity: 0});
    expect(v.getOptions(true).look['markerOpacity'], 0);
    v.setOptions({markerOpacity: 50});
    expect(v.getOptions(true).look['markerOpacity'], 50);
    v.setOptions({markerOpacity: 100});
    expect(v.getOptions(true).look['markerOpacity'], 100);
  });

  test('packCategories/yGlobalScale/axesFollowFilter/showCurrentRowLine/showMouseOverCategory combined bool round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({
      packCategories: false, yGlobalScale: true,
      axesFollowFilter: false, showCurrentRowLine: true,
      showMouseOverCategory: false,
    });
    const off = v.getOptions(true).look;
    expect(off['packCategories'], false);
    expect(off['yGlobalScale'], true);
    expect(off['axesFollowFilter'], false);
    expect(off['showCurrentRowLine'], true);
    expect(off['showMouseOverCategory'], false);
    v.setOptions({
      packCategories: true, yGlobalScale: false,
      axesFollowFilter: true, showCurrentRowLine: false,
      showMouseOverCategory: true,
    });
    const on = v.getOptions(true).look;
    expect(on['packCategories'], true);
    expect(on['yGlobalScale'], false);
    expect(on['axesFollowFilter'], true);
    expect(on['showCurrentRowLine'], false);
    expect(on['showMouseOverCategory'], true);
  });

  test('overviewColumnName + overviewAggrType combined round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.lineChart(df, {x: 'age', yColumnNames: ['height']});
    v.setOptions({overviewColumnName: 'weight', overviewAggrType: 'med'});
    const set = v.getOptions(true).look;
    expect(set['overviewColumnName'], 'weight');
    expect(set['overviewAggrType'], 'med');
    v.setOptions({overviewColumnName: '', overviewAggrType: 'avg'});
    const cleared = v.getOptions(true).look;
    expect(cleared['overviewColumnName'] === '' || cleared['overviewColumnName'] == null, true);
    expect(cleared['overviewAggrType'], 'avg');
  });
}, {owner: 'agolovko@datagrok.ai'});
