import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// JS API source: public/js-api/src/viewer.ts (DG.Viewer.boxPlot),
// public/js-api/src/interfaces/d4.d.ts:136 (IBoxPlotSettings),
// core/client/d4/lib/src/viewers/box_plot/box_plot_look.dart (BoxPlotLook).
// Box-plot-only JSON-shape coverage for props the existing
// src/ai/viewers/box-plot-js-api.ts (factory + statisticsFormat) and
// reported-issues regression tests don't pin: plotStyle (box/violin),
// violin bins int, the per-statistic boolean togglers
// (showMin/showMax/showAvg/showMed/showStdev/showQ1/showQ3),
// binColorAggrType + choices, axisType (AxisType linear/logarithmic),
// markerOpacity 0..100 boundary, and the showInside/Outside/PValue/
// EmptyCategories combined bool envelope. All assertions read state via
// getOptions(true).look — no first-paint geometry.
category('AI: Viewers: BoxPlot extras', () => {
  // no negative case: BoxPlot setOptions has no defined failure mode for
  // these look props — invalid values are coerced or ignored, not thrown.

  test('plotStyle box/violin round-trip + getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({plotStyle: 'violin'});
    expect(v.getOptions(true).look['plotStyle'], 'violin');
    v.setOptions({plotStyle: 'box'});
    expect(v.getOptions(true).look['plotStyle'], 'box');
    const props = v.getProperties();
    var styleProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'plotStyle') { styleProp = p; break; }
    expect(styleProp != null, true);
    expect(Array.isArray(styleProp!.choices), true);
    expect(styleProp!.choices.indexOf('box') >= 0, true);
    expect(styleProp!.choices.indexOf('violin') >= 0, true);
  });

  test('violin bins int round-trip within 20-1000', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({plotStyle: 'violin', bins: 100});
    expect(v.getOptions(true).look['bins'], 100);
    v.setOptions({bins: 20});
    expect(v.getOptions(true).look['bins'], 20);
    v.setOptions({bins: 1000});
    expect(v.getOptions(true).look['bins'], 1000);
  });

  test('per-statistic toggler bools combined round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({
      showStatistics: true,
      showMin: false, showMax: false, showAvg: false,
      showMed: false, showStdev: true, showQ1: true, showQ3: true,
    });
    const off = v.getOptions(true).look;
    expect(off['showMin'], false);
    expect(off['showMax'], false);
    expect(off['showAvg'], false);
    expect(off['showMed'], false);
    expect(off['showStdev'], true);
    expect(off['showQ1'], true);
    expect(off['showQ3'], true);
    v.setOptions({
      showMin: true, showMax: true, showAvg: true,
      showMed: true, showStdev: false, showQ1: false, showQ3: false,
    });
    const on = v.getOptions(true).look;
    expect(on['showMin'], true);
    expect(on['showMax'], true);
    expect(on['showAvg'], true);
    expect(on['showMed'], true);
    expect(on['showStdev'], false);
    expect(on['showQ1'], false);
    expect(on['showQ3'], false);
  });

  test('binColorAggrType round-trip + getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({binColorColumnName: 'height', binColorAggrType: 'med'});
    const look = v.getOptions(true).look;
    expect(look['binColorColumnName'], 'height');
    expect(look['binColorAggrType'], 'med');
    const props = v.getProperties();
    var aggrProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'binColorAggrType') { aggrProp = p; break; }
    expect(aggrProp != null, true);
    expect(Array.isArray(aggrProp!.choices), true);
    expect(aggrProp!.choices.length > 0, true);
    expect(aggrProp!.choices.indexOf('avg') >= 0, true);
    expect(aggrProp!.choices.indexOf('med') >= 0, true);
  });

  test('axisType linear/logarithmic round-trip + getProperties choices', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({axisType: 'logarithmic'});
    expect(v.getOptions(true).look['axisType'], 'logarithmic');
    v.setOptions({axisType: 'linear'});
    expect(v.getOptions(true).look['axisType'], 'linear');
    const props = v.getProperties();
    var axisProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'axisType') { axisProp = p; break; }
    expect(axisProp != null, true);
    expect(Array.isArray(axisProp!.choices), true);
    expect(axisProp!.choices.indexOf('linear') >= 0, true);
    expect(axisProp!.choices.indexOf('logarithmic') >= 0, true);
  });

  test('markerOpacity 0..100 boundary round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({markerOpacity: 0});
    expect(v.getOptions(true).look['markerOpacity'], 0);
    v.setOptions({markerOpacity: 50});
    expect(v.getOptions(true).look['markerOpacity'], 50);
    v.setOptions({markerOpacity: 100});
    expect(v.getOptions(true).look['markerOpacity'], 100);
  });

  test('showInsideValues/showOutsideValues/showPValue/showEmptyCategories combined bool round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.boxPlot(df, {value: 'age', category1: 'race'});
    v.setOptions({
      showInsideValues: false, showOutsideValues: false,
      showPValue: false, showEmptyCategories: false,
    });
    const off = v.getOptions(true).look;
    expect(off['showInsideValues'], false);
    expect(off['showOutsideValues'], false);
    expect(off['showPValue'], false);
    expect(off['showEmptyCategories'], false);
    v.setOptions({
      showInsideValues: true, showOutsideValues: true,
      showPValue: true, showEmptyCategories: true,
    });
    const on = v.getOptions(true).look;
    expect(on['showInsideValues'], true);
    expect(on['showOutsideValues'], true);
    expect(on['showPValue'], true);
    expect(on['showEmptyCategories'], true);
  });
}, {owner: 'agolovko@datagrok.ai'});
