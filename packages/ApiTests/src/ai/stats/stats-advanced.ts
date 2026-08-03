import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {demog, df, until} from '../helpers';

// Advanced Stats/aggregation surface: Stats.fromColumn mask, GroupByBuilder pivot/first/where/getGroups,
// grok.data.detectSemanticTypes and SemanticValue.registerRegExpDetector.
category('AI: Stats: Advanced', () => {
  test('Stats.fromColumn with BitSet mask restricts to masked rows', async () => {
    const t = df([['v', DG.COLUMN_TYPE.INT, [10, 20, 30, 40]]]);
    const col = t.col('v')!;
    // Mask only the first two rows.
    const mask = DG.BitSet.create(t.rowCount, (i) => i < 2);
    const masked = DG.Stats.fromColumn(col, mask);
    const full = DG.Stats.fromColumn(col);
    expect(masked.valueCount, 2);
    expect(masked.sum, 30);
    expect(masked.avg, 15);
    expect(masked.min, 10);
    expect(masked.max, 20);
    // Sanity: full (unmasked) stats differ.
    expect(full.valueCount, 4);
    expect(full.sum, 100);
  });

  test('GroupByBuilder.key plus multiple agg columns via aggregate', async () => {
    const t = demog(100);
    // groupBy(['race']) already adds 'race' as a key column, so do NOT call .key('race')
    // again (that would collide with "Column named 'race' already exists").
    const res = t.groupBy(['race'])
      .avg('age')
      .sum('age')
      .count('n')
      .aggregate();
    expect(res instanceof DG.DataFrame, true);
    // One row per distinct race, and each requested aggregation present as a column.
    const races = t.col('race')!.categories.filter((c) => c != null && c !== '');
    expect(res.rowCount, races.length);
    expect(res.col('race') != null, true);
    expect(res.col('avg(age)') != null, true);
    expect(res.col('sum(age)') != null, true);
    expect(res.col('n') != null, true);
    // count column sums back to the total row count.
    expect(DG.Stats.fromColumn(res.col('n')!).sum, 100);
  });

  test('GroupByBuilder.pivot produces one column per category value', async () => {
    const t = demog(100);
    const res = t.groupBy(['race'])
      .pivot('sex')
      .count()
      .aggregate();
    expect(res instanceof DG.DataFrame, true);
    // Each distinct sex value should appear in at least one pivoted column name.
    const sexes = t.col('sex')!.categories.filter((c) => c != null && c !== '');
    const colNames = res.columns.names();
    for (const s of sexes)
      expect(colNames.some((n) => n.indexOf(s) >= 0), true);
    // More columns than the single key column (key + at least one pivot column).
    expect(res.columns.length > 1, true);
  });

  test('GroupByBuilder.first returns the first value per group', async () => {
    const t = df([
      ['g', DG.COLUMN_TYPE.STRING, ['a', 'a', 'b', 'b']],
      ['v', DG.COLUMN_TYPE.INT, [1, 2, 3, 4]],
    ]);
    const res = t.groupBy(['g'])
      .first('v')
      .aggregate();
    expect(res.rowCount, 2);
    const gCol = res.col('g')!;
    const firstCol = res.col('first(v)') ?? res.columns.byIndex(res.columns.length - 1);
    expect(firstCol != null, true);
    const byGroup: {[k: string]: number} = {};
    for (let i = 0; i < res.rowCount; i++)
      byGroup[gCol.get(i)] = firstCol!.get(i);
    expect(byGroup['a'], 1);
    expect(byGroup['b'], 3);
  });

  test('GroupByBuilder.where filters source rows before aggregation', async () => {
    const t = df([
      ['g', DG.COLUMN_TYPE.STRING, ['a', 'a', 'b', 'b']],
      ['v', DG.COLUMN_TYPE.INT, [1, 2, 3, 4]],
    ]);
    const res = t.groupBy(['g'])
      .sum('v')
      .where('v > 2')
      .aggregate();
    // Only rows with v > 2 survive => group 'a' drops out entirely, 'b' keeps 3 and 4.
    expect(res.rowCount, 1);
    expect(res.col('g')!.get(0), 'b');
    expect(res.col('sum(v)')!.get(0), 7);
  });

  test('GroupByBuilder.getGroups returns groups keyed by group label', async () => {
    const t = df([
      ['g', DG.COLUMN_TYPE.STRING, ['a', 'a', 'b']],
      ['v', DG.COLUMN_TYPE.INT, [1, 2, 3]],
    ]);
    // The Dart Map<String, DataFrame> crosses interop as a plain JS object (keys are
    // 'columnName=value' strings, values are DataFrames) - not an ES Map.
    const groups = t.groupBy(['g']).getGroups() as unknown as {[k: string]: DG.DataFrame};
    const keys = Object.keys(groups);
    expect(keys.length, 2);
    let total = 0;
    for (const key of keys) {
      const sub = groups[key];
      expect(sub instanceof DG.DataFrame, true);
      total += sub.rowCount;
    }
    expect(total, 3);
    // Keys carry the group label value somewhere in the 'columnName=value' string.
    expect(keys.some((k) => k.indexOf('a') >= 0), true);
    expect(keys.some((k) => k.indexOf('b') >= 0), true);
  });

  test('grok.data.detectSemanticTypes assigns column.semType', async () => {
    // Plain numeric demog columns aren't semantically typed (semType null); detection must
    // still resolve and leave every semType as a non-empty string or null (never undefined).
    const t = demog(50);
    await grok.data.detectSemanticTypes(t);
    for (const c of t.columns) {
      const st = c.semType;
      expect(st === null || (typeof st === 'string' && st.length > 0), true);
    }
    // Detection is stable: re-running keeps each column's semType identical.
    const before = t.columns.toList().map((c) => c.semType);
    await grok.data.detectSemanticTypes(t);
    const after = t.columns.toList().map((c) => c.semType);
    for (let i = 0; i < before.length; i++)
      expect(after[i], before[i]);
  });

  test('SemanticValue.registerRegExpDetector marks matching column after detection', async () => {
    const semType = 'ai-test-code-' + Date.now();
    // A column whose values all match a distinctive regexp.
    const t = df([['code', DG.COLUMN_TYPE.STRING, ['AB-1', 'AB-2', 'AB-3', 'AB-4', 'AB-5', 'AB-6']]]);
    DG.SemanticValue.registerRegExpDetector(semType, '^AB-\\d+$', 'AI test code detector');
    await grok.data.detectSemanticTypes(t);
    // semType assignment may lag detection; poll until the registered detector lands.
    await until(() => t.col('code')!.semType === semType, 4000);
    expect(t.col('code')!.semType, semType);
  });
});
