import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {ScalarTarget, ColumnTarget, matchColumnTargets} from '../components/RunComparison/comparison-core';
import {
  ComparisonEntry, RUN_COLUMN, entryFromDataFrame,
  buildScalarComparison, buildColumnComparison, buildMultiColumnComparison,
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

category('RunComparison: multi and split dataframes', () => {
  const makeDf = (name: string, cols: DG.Column[]) => {
    const df = DG.DataFrame.fromColumns(cols);
    df.name = name;
    return df;
  };
  const floatCol = (name: string, values: number[]) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
  const intCol = (name: string, values: number[]) => DG.Column.fromList(DG.COLUMN_TYPE.INT, name, values);
  const strCol = (name: string, values: (string | null)[]) =>
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, values);
  const indexesFor = (entries: ComparisonEntry[], column: string) =>
    new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, column]])]));

  test('multi builder pads runs missing a target with nulls', async () => {
    const e1 = entryFromDataFrame(makeDf('r1', [
      intCol('time', [1, 2]), floatCol('a', [10, 20]), floatCol('b', [1, 2]),
    ]));
    const e2 = entryFromDataFrame(makeDf('r2', [
      intCol('time', [1, 2, 3]), floatCol('a', [30, 40, 50]),
    ]));
    const entries = [e1, e2];
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, 'time'))
      .filter((t): t is ColumnTarget => t.kind === 'column');
    const targetA = targets.find((t) => t.displayName === 'a')!;
    // target 'b' exists only in r1; matched targets require coverage >= 2, so bind it manually
    const targetB: ColumnTarget = {
      ...targetA,
      key: 'column:b',
      displayName: 'b',
      bindings: [{...targetA.bindings.find((x) => x.entryId === e1.id)!, columnName: 'b'}],
    };
    const result = buildMultiColumnComparison([targetA, targetB], entries)!;
    expect(result.chartDf.rowCount, 5);
    const aValues = result.chartDf.getCol('a').toList();
    const bValues = result.chartDf.getCol('b').toList();
    expect(aValues.join(','), '10,20,30,40,50');
    expect(bValues.slice(0, 2).join(','), '1,2');
    expect(bValues.slice(2).every((v: any) => v == null || isNaN(v)), true);
  });

  test('split column named Run is renamed and nulls become empty strings', async () => {
    const cols = () => [
      intCol('time', [1, 1]), floatCol('v', [10, 20]), strCol('Run', ['x', null]),
    ];
    const entries = [entryFromDataFrame(makeDf('r1', cols())), entryFromDataFrame(makeDf('r2', cols()))];
    const splits = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'Run']])]));
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, 'time'), splits);
    const result = buildColumnComparison(target as ColumnTarget, entries)!;
    expect(result.splitColumnName, 'Run (split)');
    const splitValues = result.chartDf.getCol('Run (split)').toList();
    expect(splitValues.join('|'), 'x||x|');
    expect(result.chartDf.getCol(RUN_COLUMN).toList().every((v: string) => v === 'r1' || v === 'r2'), true);
  });

  test('index kind: datetime, float, and mixed degrading to key', async () => {
    const dtCol = () => DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'ts', [dayjs('2026-01-01'), dayjs('2026-01-02')]);
    const dtEntries = [
      entryFromDataFrame(makeDf('r1', [dtCol(), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [dtCol(), floatCol('v', [3, 4])])),
    ];
    const [dtTarget] = matchColumnTargets(dtEntries.map((e) => e.nodes), indexesFor(dtEntries, 'ts'));
    const dtResult = buildColumnComparison(dtTarget as ColumnTarget, dtEntries)!;
    expect(dtResult.isKeyIndex, false);
    expect(dtResult.chartDf.getCol('ts').type, DG.COLUMN_TYPE.DATE_TIME);

    const floatEntries = [
      entryFromDataFrame(makeDf('r1', [floatCol('x', [1, 2]), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('x', [1, 2]), floatCol('v', [3, 4])])),
    ];
    const [floatTarget] = matchColumnTargets(floatEntries.map((e) => e.nodes), indexesFor(floatEntries, 'x'));
    const floatResult = buildColumnComparison(floatTarget as ColumnTarget, floatEntries)!;
    expect(floatResult.isKeyIndex, false);
    expect(floatResult.chartDf.getCol('x').type, DG.COLUMN_TYPE.FLOAT);

    const mixedEntries = [
      entryFromDataFrame(makeDf('r1', [dtCol().clone(), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('ts', [1, 2]), floatCol('v', [3, 4])])),
    ];
    const [mixedTarget] = matchColumnTargets(mixedEntries.map((e) => e.nodes), indexesFor(mixedEntries, 'ts'));
    const mixedResult = buildColumnComparison(mixedTarget as ColumnTarget, mixedEntries)!;
    expect(mixedResult.isKeyIndex, true);
    expect(mixedResult.chartDf.getCol('ts').type, DG.COLUMN_TYPE.STRING);
  });
});
