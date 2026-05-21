import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, expectNoThrow, findProp} from '../helpers';

// Regression coverage for GROK-17498: Pie chart `includeNulls` flag.
// Pins the property key on `IPieChartSettings`, the round-trip through
// `props[...]` and `getOptions(true).look`, that the viewer instance
// survives `setOptions({categoryColumnName: <other>})` while keeping both
// `includeNulls` and the new `categoryColumnName` coherent in the look
// bag, and (best-effort) that the runtime descriptor list carries
// `includeNulls` with a non-empty description.
category('AI: GROK-17498: Pie chart includeNulls property', () => {
  // Build a small df with an explicit DG missing value in the category column;
  // fall back to demog(50) (which has nulls in some columns) if the client
  // coerces JS `null` to the literal string `'null'`.
  function makeDf(): DG.DataFrame {
    const cat = DG.Column.fromList('string', 'cat', ['A', null as any, 'B', null as any, 'A']);
    if (cat.stats.missingValueCount > 0)
      return DG.DataFrame.fromColumns([cat, DG.Column.fromList('int', 'v', [1, 2, 3, 4, 5])]);
    return demog();
  }

  test('includeNulls round-trips through props and look across toggle steps', async () => {
    const df = makeDf();
    const catName = df.col('cat') != null ? 'cat' : 'race';
    expect(df.col(catName)!.stats.missingValueCount > 0, true);
    expectBoolToggle(DG.Viewer.pieChart(df, {category: catName}), 'includeNulls');
  });

  test('includeNulls survives setOptions({categoryColumnName}) switch', async () => {
    const v = DG.Viewer.pieChart(demog(), {category: 'race'});
    v.props['includeNulls'] = true;
    expectNoThrow(() => v.setOptions({categoryColumnName: 'sex'}));
    expectLook(v, {categoryColumnName: 'sex', includeNulls: true});
    v.props['includeNulls'] = false;
    expectLook(v, {categoryColumnName: 'sex', includeNulls: false});
  });

  test('includeNulls property is discoverable via getProperties (best-effort)', async () => {
    const v = DG.Viewer.pieChart(demog(20), {category: 'race'});
    const p = findProp(v, 'includeNulls');
    if (p != null && typeof p.description === 'string')
      expect(p.description.trim().length > 0, true);
    else
      expectBoolToggle(v, 'includeNulls', [true]);
  });
});
