import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

category('AI: Viewers: Histogram extras', () => {
  test('clear splitColumnName via setOptions empty string', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    v.setOptions({splitColumnName: 'race'});
    expect(v.props['splitColumnName'], 'race');
    expect(v.getOptions(true).look['splitColumnName'], 'race');
    v.setOptions({splitColumnName: ''});
    const cleared = v.getOptions(true).look['splitColumnName'];
    const propCleared = v.props['splitColumnName'];
    expect(cleared === '' || cleared == null, true);
    expect(propCleared === '' || propCleared == null, true);
    expect(cleared, propCleared);
  });

  test('filter formula round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const formula = '${age} > 20';
    v.setOptions({filter: formula});
    expect(v.props['filter'], formula);
    expect(v.getOptions(true).look['filter'], formula);
  });

  test('linearColorScheme number-array round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const scheme = [0xFF112233, 0xFFAABBCC];
    v.setOptions({linearColorScheme: scheme});
    const look = v.getOptions(true).look;
    expect(Array.isArray(look['linearColorScheme']), true);
    expect((look['linearColorScheme'] as number[]).length, 2);
    expectArray(look['linearColorScheme'] as number[], scheme);
  });

  test('getProperties columnTypeFilter introspection', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const props = v.getProperties();
    const byName: {[k: string]: DG.Property} = {};
    for (var p of props)
      byName[p.name] = p;
    expect(byName['valueColumnName'] != null, true);
    expect(byName['splitColumnName'] != null, true);
    expect(byName['colorColumnName'] != null, true);
    expect(byName['valueColumnName'].columnTypeFilter as string, 'numerical');
    expect(byName['splitColumnName'].columnTypeFilter as string, 'categorical');
    expect(byName['colorColumnName'].columnTypeFilter as string, 'numerical');
  });

  test('getProperties choices introspection on colorAggrType', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const props = v.getProperties();
    var aggrProp: DG.Property | null = null;
    for (var p of props)
      if (p.name === 'colorAggrType') { aggrProp = p; break; }
    expect(aggrProp != null, true);
    expect(Array.isArray(aggrProp!.choices), true);
    expect(aggrProp!.choices.length > 0, true);
  });

  test('aggTooltipColumns string round-trip', async () => {
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.histogram(df, {value: 'age'});
    const tooltip = 'avg(age), count(race)';
    v.setOptions({aggTooltipColumns: tooltip});
    expect(v.props['aggTooltipColumns'], tooltip);
    expect(v.getOptions(true).look['aggTooltipColumns'], tooltip);
  });
});
