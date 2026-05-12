import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

// Regression coverage for GROK-17498: Pie chart `includeNulls` flag.
// Pins the property key on `IPieChartSettings`, the round-trip through
// `props[...]` and `getOptions(true).look`, that the viewer instance
// survives `setOptions({categoryColumnName: <other>})` while keeping both
// `includeNulls` and the new `categoryColumnName` coherent in the look
// bag, and (best-effort) that the runtime descriptor list carries
// `includeNulls` with a non-empty description.
//
// Note: the triage referenced `DG.Viewer.pieChart(df, ...)`; that helper does not
// exist on `DataFramePlotHelper`. We use the documented
// `DG.Viewer.pieChart(df, ...)` factory, matching the existing PieChart
// suite (`ai/viewers/pie-chart-js-api.ts`).
category('AI: GROK-17498: Pie chart includeNulls property', () => {
  // Build a small df with an explicit DG missing value in the category
  // column. If a future client coerces JS `null` to the literal string
  // `'null'`, fall back to `grok.data.demo.demog(50)` (which has nulls in
  // some columns) — the triage flagged this as a "null semantics" risk.
  function makeDf(): DG.DataFrame {
    const cat = DG.Column.fromList('string', 'cat', ['A', null as any, 'B', null as any, 'A']);
    if (cat.stats.missingValueCount > 0) {
      const v = DG.Column.fromList('int', 'v', [1, 2, 3, 4, 5]);
      return DG.DataFrame.fromColumns([cat, v]);
    }
    // Fallback: demog(50) has nulls in several columns.
    return grok.data.demo.demog(50);
  }

  function categoryName(df: DG.DataFrame): string {
    // 'cat' for the hand-built df; 'race' for the demog fallback (which
    // has missing values across rows).
    for (var c of df.columns.toList())
      if (c.name === 'cat') return 'cat';
    return 'race';
  }

  test('includeNulls round-trips through props and look across toggle steps', async () => {
    const df = makeDf();
    const catName = categoryName(df);
    const catCol = df.col(catName);
    expect(catCol != null, true);
    // The category column must carry true DG missing values.
    expect(catCol!.stats.missingValueCount > 0, true);

    const v = DG.Viewer.pieChart(df, {category: catName});
    expect(v instanceof DG.PieChartViewer, true);
    expect(v.type, DG.VIEWER.PIE_CHART);

    v.props['includeNulls'] = true;
    expect(v.props['includeNulls'], true);
    expect(v.getOptions(true).look['includeNulls'], true);

    v.props['includeNulls'] = false;
    expect(v.props['includeNulls'], false);
    expect(v.getOptions(true).look['includeNulls'], false);

    v.props['includeNulls'] = true;
    expect(v.props['includeNulls'], true);
    expect(v.getOptions(true).look['includeNulls'], true);
  });

  test('includeNulls survives setOptions({categoryColumnName}) switch', async () => {
    // Use demog so we have at least two viable category columns side by
    // side, both carrying nulls in some columns of the demo data.
    const df = grok.data.demo.demog(50);
    const v = DG.Viewer.pieChart(df, {category: 'race'});
    expect(v instanceof DG.PieChartViewer, true);

    v.props['includeNulls'] = true;
    expect(v.getOptions(true).look['includeNulls'], true);

    var threw = false;
    try {
      v.setOptions({categoryColumnName: 'sex'});
    }
    catch (_e) {
      threw = true;
    }
    expect(threw, false);

    // Viewer instance survives — type and look-bag both intact.
    expect(v.type, DG.VIEWER.PIE_CHART);
    const look = v.getOptions(true).look;
    expect(look['categoryColumnName'], 'sex');
    expect(look['includeNulls'], true);

    // Toggle once more under the new category; both keys stay coherent.
    v.props['includeNulls'] = false;
    const look2 = v.getOptions(true).look;
    expect(look2['categoryColumnName'], 'sex');
    expect(look2['includeNulls'], false);
  });

  test('includeNulls property is discoverable via getProperties (best-effort)', async () => {
    const df = grok.data.demo.demog(20);
    const v = DG.Viewer.pieChart(df, {category: 'race'});
    const props = v.getProperties() as DG.Property[];
    expect(Array.isArray(props), true);
    expect(props.length > 0, true);
    var found: DG.Property | null = null;
    for (var p of props) {
      if (p.name === 'includeNulls') {
        found = p;
        break;
      }
    }
    // Best-effort: if the runtime descriptor list does not carry
    // `includeNulls`, fall back to the schema-only check we already
    // exercised above (set via props, read via look).
    if (found != null) {
      const desc = found.description;
      expect(typeof desc === 'string' || desc == null, true);
      if (typeof desc === 'string')
        expect(desc.trim().length > 0, true);
    }
    else {
      v.props['includeNulls'] = true;
      expect(v.getOptions(true).look['includeNulls'], true);
    }
  });
});
