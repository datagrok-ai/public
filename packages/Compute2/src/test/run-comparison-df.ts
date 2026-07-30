import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {ScalarTarget, matchColumnTargets} from '../components/RunComparison/comparison-core';
import {
  ComparisonEntry, RUN_COLUMN, entryFromDataFrame, buildScalarComparison, buildColumnComparison,
} from '../components/RunComparison/data-extraction';

// Regression: run/index columns must be created with an enforced string type.
// Column.fromStrings infers the type from values, so numeric-looking run names
// used to produce a numeric column and break chart legends (getUsedCategories).

function scalarEntry(id: string, name: string, value: number): ComparisonEntry {
  return {
    id,
    name,
    sourceKind: 'function',
    modelName: '',
    nodes: {
      entryId: id,
      entryName: name,
      scalars: [{path: 'x', name: 'x', valueType: 'double', value}],
      tables: [],
    },
    dataFrames: new Map(),
  };
}

category('RunComparison: comparison dataframes', () => {
  test('run column stays string for numeric-looking run names', async () => {
    const entries = [scalarEntry('a', '1', 10), scalarEntry('b', '2', 20)];
    const target: ScalarTarget = {
      kind: 'scalar',
      key: 'scalar:x',
      displayName: 'x',
      confidence: 'exact',
      unitsWarning: false,
      coverage: 2,
      total: 2,
      bindings: [
        {entryId: 'a', path: 'x', name: 'x', value: 10},
        {entryId: 'b', path: 'x', name: 'x', value: 20},
      ],
    };
    const {chartDf} = buildScalarComparison(target, entries);
    expect(chartDf.getCol(RUN_COLUMN).type, DG.COLUMN_TYPE.STRING);
    expect(chartDf.getCol('Path').type, DG.COLUMN_TYPE.STRING);
  });

  test('keyed index and run columns stay string in column comparison', async () => {
    const makeDf = (name: string) => {
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'key', ['1', '2', '3']),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20, 30]),
      ]);
      df.name = name;
      return df;
    };
    // numeric-looking entry names (raw table names)
    const entries = [entryFromDataFrame(makeDf('1')), entryFromDataFrame(makeDf('2'))];
    const indexes = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'key']])]));
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    const result = buildColumnComparison(target, entries)!;
    expect(result.isKeyIndex, true);
    expect(result.chartDf.getCol('key').type, DG.COLUMN_TYPE.STRING);
    expect(result.chartDf.getCol(RUN_COLUMN).type, DG.COLUMN_TYPE.STRING);
  });
});
