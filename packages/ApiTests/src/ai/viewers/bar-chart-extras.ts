import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: BarChart extras', () => {
  test('color friendly-key alias on factory', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race', color: 'height'});
    expect(v.type, DG.VIEWER.BAR_CHART);
    expect(v instanceof DG.BarChartViewer, true);
    expect(v.props['colorColumnName'], 'height');
    const look = v.getOptions(true).look;
    expect(look['colorColumnName'], 'height');
  });

  test('split and stack appear together in look and props', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race', stack: 'sex'});
    expect(v.props['splitColumnName'], 'race');
    expect(v.props['stackColumnName'], 'sex');
    const look = v.getOptions(true).look;
    expect(look['splitColumnName'], 'race');
    expect(look['stackColumnName'], 'sex');
  });

  test('barSortType and barSortOrder round-trip together', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'});
    v.setOptions({barSortType: 'by value', barSortOrder: 'asc'});
    expect(v.props['barSortType'], 'by value');
    expect(v.props['barSortOrder'], 'asc');
    const look = v.getOptions(true).look;
    expect(look['barSortType'], 'by value');
    expect(look['barSortOrder'], 'asc');
  });

  test('getProperties choices introspection on barSortType', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'});
    const props = v.getProperties();
    var sortProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'barSortType') { sortProp = p; break; }
    expect(sortProp != null, true);
    expect(Array.isArray(sortProp!.choices), true);
    expect(sortProp!.choices.length > 0, true);
    expect(sortProp!.choices.indexOf('by category') >= 0, true);
    expect(sortProp!.choices.indexOf('by value') >= 0, true);
  });

  test('includeNulls boolean round-trip via setOptions', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'});
    v.setOptions({includeNulls: false});
    expect(v.props['includeNulls'], false);
    expect(v.getOptions(true).look['includeNulls'], false);
    v.setOptions({includeNulls: true});
    expect(v.props['includeNulls'], true);
    expect(v.getOptions(true).look['includeNulls'], true);
  });

  test('onClick RowGroupAction enum round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.barChart(df, {value: 'age', split: 'race'});
    v.setOptions({onClick: 'Filter'});
    expect(v.props['onClick'], 'Filter');
    expect(v.getOptions(true).look['onClick'], 'Filter');
    v.setOptions({onClick: 'Select'});
    expect(v.props['onClick'], 'Select');
    expect(v.getOptions(true).look['onClick'], 'Select');
  });
});
