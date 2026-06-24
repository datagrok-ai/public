import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {
  demog, expectChoices, expectPropAndLook, expectRoundTripPropAndLook, findProp, subscribeAll, withAttachedViewer,
} from '../helpers';

// PieChart JS API: factory, category alias, look round-trip, events, addViewer.
category('AI: Viewers: PieChart JS API', () => {
  test('factory typed', async () => {
    const v = DG.Viewer.pieChart(demog(100), {category: 'race'});
    expect(v instanceof DG.PieChartViewer, true);
    expectPropAndLook(v, {categoryColumnName: 'race'});
  });

  test('pie-specific look round-trip via setOptions', async () => {
    const v = DG.Viewer.pieChart(demog(), {category: 'race'});
    expectRoundTripPropAndLook(v, {
      pieSortType: 'by category', pieSortOrder: 'desc', includeNulls: false, labelPosition: 'Outside',
    });
  });

  test('segment angle and length columns', async () => {
    const v = DG.Viewer.pieChart(demog(), {category: 'race'});
    expectRoundTripPropAndLook(v, {
      segmentAngleColumnName: 'age', segmentAngleAggrType: 'avg',
      segmentLengthColumnName: 'weight', segmentLengthAggrType: 'max',
    });
    v.setOptions({segmentAngleColumnName: '', segmentLengthColumnName: ''});
    expect(v.props['segmentAngleColumnName'], '');
    expect(v.props['segmentLengthColumnName'], '');
  });

  test('getProperties exposes pie-specific properties', async () => {
    const v = DG.Viewer.pieChart(demog(20), {category: 'race'});
    for (const name of ['categoryColumnName', 'pieSortType', 'pieSortOrder',
      'segmentAngleColumnName', 'segmentLengthColumnName', 'includeNulls'])
      expect(findProp(v, name) != null, true);
    expectChoices(v, 'pieSortType', ['by value', 'by category']);
  });

  test('onSegmentClicked is an rxjs Observable', async () => {
    const v = DG.Viewer.pieChart(demog(20), {category: 'race'}) as DG.PieChartViewer;
    subscribeAll([v.onSegmentClicked])();
  });

  test('view.addViewer attaches a typed PieChartViewer', async () => {
    await withAttachedViewer<DG.PieChartViewer>(demog(20), DG.VIEWER.PIE_CHART, {category: 'race'},
      (v) => expect(v instanceof DG.PieChartViewer, true));
  });
});
