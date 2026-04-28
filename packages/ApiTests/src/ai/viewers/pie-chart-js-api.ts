import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: PieChart JS API', () => {
  test('factory typed', async () => {
    const df = grok.data.demo.demog(100);
    const v = DG.Viewer.pieChart(df, {category: 'race'});
    expect(v instanceof DG.PieChartViewer, true);
    expect(v.type, DG.VIEWER.PIE_CHART);
    expect(v.props['categoryColumnName'], 'race');
    const look = v.getOptions(true).look;
    expect(look['categoryColumnName'], 'race');
  });

  test('pie-specific look round-trip via setOptions', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pieChart(df, {category: 'race'});
    v.setOptions({
      pieSortType: 'by category',
      pieSortOrder: 'desc',
      includeNulls: false,
      labelPosition: 'Outside',
    });
    expect(v.props['pieSortType'], 'by category');
    expect(v.props['pieSortOrder'], 'desc');
    expect(v.props['includeNulls'], false);
    expect(v.props['labelPosition'], 'Outside');
    const look = v.getOptions(true).look;
    expect(look['pieSortType'], 'by category');
    expect(look['pieSortOrder'], 'desc');
    expect(look['includeNulls'], false);
    expect(look['labelPosition'], 'Outside');
  });

  test('segment angle and length columns', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pieChart(df, {category: 'race'});
    v.setOptions({
      segmentAngleColumnName: 'age',
      segmentAngleAggrType: 'avg',
      segmentLengthColumnName: 'weight',
      segmentLengthAggrType: 'max',
    });
    expect(v.props['segmentAngleColumnName'], 'age');
    expect(v.props['segmentAngleAggrType'], 'avg');
    expect(v.props['segmentLengthColumnName'], 'weight');
    expect(v.props['segmentLengthAggrType'], 'max');
    const look = v.getOptions(true).look;
    expect(look['segmentAngleColumnName'], 'age');
    expect(look['segmentLengthColumnName'], 'weight');

    v.setOptions({segmentAngleColumnName: '', segmentLengthColumnName: ''});
    expect(v.props['segmentAngleColumnName'], '');
    expect(v.props['segmentLengthColumnName'], '');
  });

  test('getProperties exposes pie-specific properties', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pieChart(df, {category: 'race'});
    const props = v.getProperties();
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);
    const expected = [
      'categoryColumnName',
      'pieSortType',
      'pieSortOrder',
      'segmentAngleColumnName',
      'segmentLengthColumnName',
      'includeNulls',
    ];
    for (var name of expected) {
      var found = false;
      for (var p of props)
        if (p.name === name) { found = true; break; }
      expect(found, true);
    }
    var sortTypeProp: DG.Property | undefined;
    for (var p of props)
      if (p.name === 'pieSortType') { sortTypeProp = p; break; }
    expect(sortTypeProp != null, true);
    const choices = sortTypeProp!.choices;
    expect(Array.isArray(choices), true);
    expect(choices.indexOf('by value') >= 0, true);
    expect(choices.indexOf('by category') >= 0, true);
  });

  test('onSegmentClicked is an rxjs Observable', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pieChart(df, {category: 'race'}) as DG.PieChartViewer;
    expect(v instanceof DG.PieChartViewer, true);
    expect(v.onSegmentClicked != null, true);
    expect(typeof v.onSegmentClicked.subscribe, 'function');
    const subs: Subscription[] = [];
    try {
      const sub = v.onSegmentClicked.subscribe(() => {});
      expect(typeof sub.unsubscribe, 'function');
      subs.push(sub);
    }
    finally {
      for (var sub of subs)
        sub.unsubscribe();
    }
  });

  test('view.addViewer attaches a typed PieChartViewer', async () => {
    const df = grok.data.demo.demog(20);
    const tv = grok.shell.addTableView(df);
    try {
      const v = tv.addViewer(DG.VIEWER.PIE_CHART, {category: 'race'});
      expect(v.type, DG.VIEWER.PIE_CHART);
      expect(v instanceof DG.PieChartViewer, true);
      var found: DG.Viewer | undefined;
      for (var x of tv.viewers)
        if (x.type === DG.VIEWER.PIE_CHART) { found = x; break; }
      expect(found != null, true);
      expect(found instanceof DG.PieChartViewer, true);
    }
    finally {
      tv.close();
    }
  });
});
