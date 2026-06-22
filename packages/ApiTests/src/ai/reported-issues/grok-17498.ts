import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectBoolToggle, expectLook, expectNoThrow, findProp} from '../helpers';

// Regression coverage for GROK-17498: Pie chart includeNulls flag.
category('AI: GROK-17498: Pie chart includeNulls property', () => {
  // Fall back to demog if the client coerces JS null to the string 'null'.
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

  test('includeNulls property is discoverable via getProperties', async () => {
    const v = DG.Viewer.pieChart(demog(20), {category: 'race'});
    expect(findProp(v, 'includeNulls') != null, true);
  });
});
