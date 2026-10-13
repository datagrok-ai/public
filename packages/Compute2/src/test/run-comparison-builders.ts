import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {
  ScalarTarget, ColumnTarget, ComparisonEntry, RUN_COLUMN, TimeUnit,
} from '../components/RunComparison/types';
import {matchScalarTargets, matchColumnTargets, buildUnitConverter} from '../components/RunComparison/matching';
import {entryFromDataFrame} from '../components/RunComparison/entry-extraction';
import {
  buildScalarComparison, buildMultiScalarComparison, buildColumnComparison, buildComparison, columnPrimitives,
} from '../components/RunComparison/comparison-builders';

// Regression: run/index columns must be created with an enforced string type.
// Column.fromStrings infers the type from values, so numeric-looking run names
// used to produce a numeric column and break chart legends (getUsedCategories).

function scalarsEntry(id: string, name: string, scalars: {name: string, value: number}[]): ComparisonEntry {
  return {
    id,
    name,
    sourceKind: 'function',
    modelName: '',
    nodes: {
      entryId: id,
      entryName: name,
      scalars: scalars.map((s) => ({path: s.name, name: s.name, valueType: 'double', value: s.value})),
      tables: [],
    },
    dataFrames: new Map(),
  };
}

function scalarEntry(id: string, name: string, value: number): ComparisonEntry {
  return scalarsEntry(id, name, [{name: 'x', value}]);
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
    const result = buildColumnComparison([target], entries)!;
    expect(result.isKeyIndex, true);
    expect(result.chartDf.getCol('key').type, DG.COLUMN_TYPE.STRING);
    expect(result.chartDf.getCol(RUN_COLUMN).type, DG.COLUMN_TYPE.STRING);
  });

  test('missing key index cells stay missing', async () => {
    const makeDf = (name: string, keys: (string | number | null)[], type: DG.ColumnType) => {
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList(type, 'key', keys),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20]),
      ]);
      df.name = name;
      return df;
    };
    const entries = [
      entryFromDataFrame(makeDf('r1', ['a', null], DG.COLUMN_TYPE.STRING)),
      entryFromDataFrame(makeDf('r2', [1, null], DG.COLUMN_TYPE.INT)),
    ];
    const indexes = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'key']])]));
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    const key = buildColumnComparison([target], entries)!.chartDf.getCol('key');
    expectArray(key.categories.filter((c) => c !== ''), ['1', 'a']);
    expect(key.isNone(1) && key.isNone(3), true);
  });

  test('returns null when fewer than two runs participate', async () => {
    const makeDf = (name: string) => {
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.INT, 'time', [1, 2]),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20]),
      ]);
      df.name = name;
      return df;
    };
    const entries = [entryFromDataFrame(makeDf('r1')), entryFromDataFrame(makeDf('r2'))];
    const indexes = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'time']])]));
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    const single: ColumnTarget = {...(target as ColumnTarget), bindings: [target.bindings[0]]};
    expect(buildColumnComparison([single], entries), null);
  });

  test('scalar target named Run does not collide with the run column', async () => {
    const entries = [scalarsEntry('a', 'run one', [{name: 'Run', value: 10}]),
      scalarsEntry('b', 'run two', [{name: 'Run', value: 20}])];
    const [target] = matchScalarTargets(entries.map((e) => e.nodes)) as ScalarTarget[];
    const {chartDf, valueColumnName} = buildScalarComparison(target, entries);
    expect(chartDf.getCol(RUN_COLUMN).type, DG.COLUMN_TYPE.STRING);
    expectArray(chartDf.getCol(RUN_COLUMN).toList(), ['run one', 'run two']);
    // the column the chart reads its values from must be the numeric one
    expect(valueColumnName === RUN_COLUMN, false);
    expect(chartDf.getCol(valueColumnName).type, DG.COLUMN_TYPE.FLOAT);
    expectArray(chartDf.getCol(valueColumnName).toList(), [10, 20]);
  });

  test('column target named Run does not collide with the run column', async () => {
    const makeDf = (name: string) => {
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList(DG.COLUMN_TYPE.INT, 'time', [1, 2]),
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Run', [10, 20]),
      ]);
      df.name = name;
      return df;
    };
    const entries = [entryFromDataFrame(makeDf('r1')), entryFromDataFrame(makeDf('r2'))];
    const indexes = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'time']])]));
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    const result = buildColumnComparison([target], entries)!;
    expect(result.chartDf.getCol(RUN_COLUMN).type, DG.COLUMN_TYPE.STRING);
    expectArray(result.chartDf.getCol(RUN_COLUMN).toList(), ['r1', 'r1', 'r2', 'r2']);
    // the column the chart reads its values from must be the numeric one
    expect(result.valueColumnName === RUN_COLUMN, false);
    expect(result.chartDf.getCol(result.valueColumnName).type, DG.COLUMN_TYPE.FLOAT);
    expectArray(result.chartDf.getCol(result.valueColumnName).toList(), [10, 20, 10, 20]);
  });
});

category('RunComparison: multi-scalar dataframes', () => {
  const scalarTargetsFor = (entries: ComparisonEntry[]) =>
    matchScalarTargets(entries.map((e) => e.nodes)).filter((t): t is ScalarTarget => t.kind === 'scalar');

  test('one row per run, string run column, values in entry order', async () => {
    const entries = [
      scalarsEntry('a', '1', [{name: 'x', value: 10}, {name: 'y', value: 1}]),
      scalarsEntry('b', '2', [{name: 'x', value: 20}, {name: 'y', value: 2}]),
      scalarsEntry('c', '3', [{name: 'x', value: 30}, {name: 'y', value: 3}]),
    ];
    const targets = scalarTargetsFor(entries);
    const result = buildMultiScalarComparison(targets, entries)!;
    expect(result.chartDf.rowCount, 3);
    expect(result.chartDf.getCol(RUN_COLUMN).type, DG.COLUMN_TYPE.STRING);
    expect(result.chartDf.getCol(RUN_COLUMN).toList().join('|'), '1|2|3');
    expect(result.valueColumnNames.join('|'), 'x|y');
    expect(result.chartDf.getCol('x').toList().join(','), '10,20,30');
    expect(result.chartDf.getCol('y').toList().join(','), '1,2,3');
  });

  test('run missing a scalar is padded with null, unmatched run excluded', async () => {
    const entries = [
      scalarsEntry('a', 'r1', [{name: 'x', value: 10}, {name: 'y', value: 1}]),
      scalarsEntry('b', 'r2', [{name: 'x', value: 20}, {name: 'y', value: 2}]),
      scalarsEntry('c', 'r3', [{name: 'x', value: 30}]),
      scalarsEntry('d', 'r4', [{name: 'unrelated', value: 0}]),
    ];
    const targets = scalarTargetsFor(entries);
    const result = buildMultiScalarComparison(targets, entries)!;
    expect(result.chartDf.rowCount, 3);
    expect(result.chartDf.getCol(RUN_COLUMN).toList().join('|'), 'r1|r2|r3');
    const yValues = result.chartDf.getCol('y').toList();
    expect(yValues.slice(0, 2).join(','), '1,2');
    expect(yValues[2] == null || isNaN(yValues[2]), true);
  });

  test('duplicate display names get deduplicated value columns', async () => {
    const entries = [
      scalarsEntry('a', 'r1', [{name: 'x', value: 10}]),
      scalarsEntry('b', 'r2', [{name: 'x', value: 20}]),
    ];
    const [target] = scalarTargetsFor(entries);
    const twin: ScalarTarget = {...target, key: 'scalar:x:2'};
    const result = buildMultiScalarComparison([target, twin], entries)!;
    expect(result.valueColumnNames.join('|'), 'x|x (2)');
    expect(result.chartDf.getCol('x (2)').toList().join(','), '10,20');
  });

  test('returns null when fewer than two runs participate', async () => {
    const entries = [
      scalarsEntry('a', 'r1', [{name: 'x', value: 10}, {name: 'y', value: 1}]),
      scalarsEntry('b', 'r2', [{name: 'x', value: 20}, {name: 'y', value: 2}]),
    ];
    const targets = scalarTargetsFor(entries);
    const singles = targets.map((t) => ({...t, bindings: t.bindings.slice(0, 1)}));
    expect(buildMultiScalarComparison(singles, entries), null);
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
    const result = buildColumnComparison([targetA, targetB], entries)!;
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
    const result = buildColumnComparison([target as ColumnTarget], entries)!;
    expect(result.splitColumnName, 'Run (split)');
    const splitValues = result.chartDf.getCol('Run (split)').toList();
    expect(splitValues.join('|'), 'x||x|');
    expect(result.chartDf.getCol(RUN_COLUMN).toList().every((v: string) => v === 'r1' || v === 'r2'), true);
  });

  test('unsplit runs pad the split column with empty strings', async () => {
    const e1 = entryFromDataFrame(makeDf('r1', [
      intCol('time', [1, 1]), floatCol('v', [10, 20]), strCol('batch', ['A', 'B']),
    ]));
    const e2 = entryFromDataFrame(makeDf('r2', [
      intCol('time', [1, 2]), floatCol('v', [30, 40]),
    ]));
    const entries = [e1, e2];
    const splits = new Map([[e1.id, new Map([[e1.nodes.tables[0].path, 'batch']])]]);
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, 'time'), splits);
    expect((target as ColumnTarget).coverage, 2);
    const result = buildColumnComparison([target as ColumnTarget], entries)!;
    expect(result.splitColumnName, 'batch');
    expect(result.chartDf.getCol('batch').toList().join('|'), 'A|B||');
    expect(result.chartDf.getCol('v').toList().join(','), '10,20,30,40');
  });

  test('index kind: datetime, float, and mixed degrading to key', async () => {
    const dtCol = () => DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'ts', [dayjs('2026-01-01'), dayjs('2026-01-02')]);
    const dtEntries = [
      entryFromDataFrame(makeDf('r1', [dtCol(), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [dtCol(), floatCol('v', [3, 4])])),
    ];
    const [dtTarget] = matchColumnTargets(dtEntries.map((e) => e.nodes), indexesFor(dtEntries, 'ts'));
    const dtResult = buildColumnComparison([dtTarget as ColumnTarget], dtEntries)!;
    expect(dtResult.isKeyIndex, false);
    expect(dtResult.chartDf.getCol('ts').type, DG.COLUMN_TYPE.DATE_TIME);

    const floatEntries = [
      entryFromDataFrame(makeDf('r1', [floatCol('x', [1, 2]), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('x', [1, 2]), floatCol('v', [3, 4])])),
    ];
    const [floatTarget] = matchColumnTargets(floatEntries.map((e) => e.nodes), indexesFor(floatEntries, 'x'));
    const floatResult = buildColumnComparison([floatTarget as ColumnTarget], floatEntries)!;
    expect(floatResult.isKeyIndex, false);
    expect(floatResult.chartDf.getCol('x').type, DG.COLUMN_TYPE.FLOAT);

    const mixedEntries = [
      entryFromDataFrame(makeDf('r1', [dtCol().clone(), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('ts', [1, 2]), floatCol('v', [3, 4])])),
    ];
    const [mixedTarget] = matchColumnTargets(mixedEntries.map((e) => e.nodes), indexesFor(mixedEntries, 'ts'));
    const mixedResult = buildColumnComparison([mixedTarget as ColumnTarget], mixedEntries)!;
    expect(mixedResult.isKeyIndex, true);
    expect(mixedResult.chartDf.getCol('ts').type, DG.COLUMN_TYPE.STRING);
  });

  test('cross-table pick pads that run with nulls', async () => {
    const r1 = entryFromDataFrame(makeDf('r1', [
      intCol('time', [1, 2]), floatCol('a', [10, 20]), floatCol('b', [1, 2]),
    ]));
    const r2: ComparisonEntry = {
      id: 'r2',
      name: 'r2',
      sourceKind: 'function',
      modelName: '',
      nodes: {
        entryId: 'r2',
        entryName: 'r2',
        scalars: [],
        tables: [
          {path: 'p1', name: 'p1', rowCount: 3,
            columns: [{name: 'time', type: 'int'}, {name: 'a', type: 'double'}]},
          {path: 'p2', name: 'p2', rowCount: 2,
            columns: [{name: 'time', type: 'int'}, {name: 'b', type: 'double'}]},
        ],
      },
      dataFrames: new Map([
        ['p1', makeDf('p1', [intCol('time', [1, 2, 3]), floatCol('a', [30, 40, 50])])],
        ['p2', makeDf('p2', [intCol('time', [1, 2]), floatCol('b', [7, 8])])],
      ]),
    };
    const entries = [r1, r2];
    const r1Table = r1.nodes.tables[0].path;
    const bindingFor = (entryId: string, tablePath: string, columnName: string) =>
      ({entryId, tablePath, tableName: tablePath, columnName, indexColumnName: 'time'});
    const targetFor = (key: string, bindings: ReturnType<typeof bindingFor>[]): ColumnTarget => ({
      kind: 'column', key, displayName: key, confidence: 'exact', unitsWarning: false,
      coverage: 2, total: 2, candidates: [], bindings,
    });
    const anchor = targetFor('a', [bindingFor(r1.id, r1Table, 'a'), bindingFor('r2', 'p1', 'a')]);
    const co = targetFor('b', [bindingFor(r1.id, r1Table, 'b'), bindingFor('r2', 'p2', 'b')]);
    const result = buildColumnComparison([anchor, co], entries)!;
    expect(result.chartDf.rowCount, 5);
    expect(result.chartDf.getCol('a').toList().join(','), '10,20,30,40,50');
    const bValues = result.chartDf.getCol('b').toList();
    expect(bValues.slice(0, 2).join(','), '1,2');
    expect(bValues.slice(2).every((v: any) => v == null || isNaN(v)), true);
  });

  test('disabled candidates are excluded from the chart', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [intCol('time', [1, 2]), floatCol('value', [10, 20])])),
      entryFromDataFrame(makeDf('r2', [intCol('time', [1, 2]), floatCol('value', [30, 40])])),
    ];
    const indexes = indexesFor(entries, 'time');
    const [base] = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexes, undefined,
      {[base.key]: {[`${entries[1].id}|r2|value`]: false}});
    expect(buildColumnComparison([target as ColumnTarget], entries), null);
  });

  test('duplicate display names get deduplicated value columns', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [intCol('time', [1, 2]), floatCol('x', [10, 20])])),
      entryFromDataFrame(makeDf('r2', [intCol('time', [1, 2]), floatCol('x', [30, 40])])),
    ];
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, 'time'));
    const twin: ColumnTarget = {...(target as ColumnTarget), key: 'column:x:twin'};
    const result = buildColumnComparison([target as ColumnTarget, twin], entries)!;
    expect(result.valueColumnNames.join('|'), 'x|x (2)');
    expect(result.chartDf.getCol('x (2)').toList().join(','), '10,20,30,40');
  });
});

category('RunComparison: time series dataframes', () => {
  const makeDf = (name: string, cols: DG.Column[]) => {
    const df = DG.DataFrame.fromColumns(cols);
    df.name = name;
    return df;
  };
  const floatCol = (name: string, values: (number | null)[]) =>
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
  const dtCol = (name: string, values: string[]) =>
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, name, values.map((v) => dayjs(v)));
  const indexesFor = (entries: ComparisonEntry[], column: string) =>
    new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, column]])]));
  const targetFor = (entries: ComparisonEntry[], column: string) =>
    matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, column))
      .find((t): t is ColumnTarget => t.kind === 'column')!;
  const tsFor = (entries: ComparisonEntry[], units: (TimeUnit | undefined)[]) =>
    new Map(entries.map((e, i) => [e.id,
      new Map([[e.nodes.tables[0].path, units[i] === undefined ?
        {mode: 'timeseries' as const} : {mode: 'timeseries' as const, units: units[i]}]])]));
  // Array.from: toList() returns sparse arrays for null cells, and map() would skip the holes
  const numbers = (values: any[]) => Array.from(values, (v) => v == null || isNaN(v) ? null : v);

  test('aligned numeric indexes with heterogeneous units share one axis', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [10, 20]), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [1, 2]), floatCol('v', [3, 4])])),
    ];
    const result = buildColumnComparison([targetFor(entries, 'time')], entries,
      tsFor(entries, ['s', 'min']))!;
    expect(result.isKeyIndex, false);
    expect(result.timeSeriesUnit, 's');
    expect(result.indexColumnName, 'time (s)');
    expect(result.chartDf.getCol('time (s)').type, DG.COLUMN_TYPE.FLOAT);
    expect(numbers(result.chartDf.getCol('time (s)').toList()).join(','), '0,10,0,60');
  });

  test('bigint indexes convert through the elapsed path', async () => {
    const bigintCol = (name: string, values: bigint[]) =>
      DG.Column.fromBigInt64Array(name, new BigInt64Array(values));
    const entries = [
      entryFromDataFrame(makeDf('r1', [bigintCol('time', [BigInt(10), BigInt(20)]), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [bigintCol('time', [BigInt(5), BigInt(15)]), floatCol('v', [3, 4])])),
    ];
    const result = buildColumnComparison([targetFor(entries, 'time')], entries,
      tsFor(entries, ['s', 's']))!;
    expect(result.timeSeriesUnit, 's');
    expect(numbers(result.chartDf.getCol('time (s)').toList()).join(','), '0,10,0,10');
  });

  test('unsorted index aligns to its minimum', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [30, 10, 20]), floatCol('v', [1, 2, 3])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [5, 15]), floatCol('v', [4, 5])])),
    ];
    const result = buildColumnComparison([targetFor(entries, 'time')], entries,
      tsFor(entries, ['s', 's']))!;
    expect(numbers(result.chartDf.getCol('time (s)').toList()).join(','), '20,0,10,0,10');
  });

  test('mixed datetime and float tables chart as an aligned line', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [dtCol('ts', ['2026-01-01', '2026-01-02']), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('ts', [0, 60]), floatCol('v', [3, 4])])),
    ];
    const aligned = buildColumnComparison([targetFor(entries, 'ts')], entries,
      tsFor(entries, [undefined, 'min']))!;
    expect(aligned.isKeyIndex, false);
    expect(aligned.timeSeriesUnit, 'min');
    expect(numbers(aligned.chartDf.getCol('ts (min)').toList()).join(','), '0,1440,0,60');
  });

  test('datetime-only comparison picks an auto unit from the range', async () => {
    const cols = () => [dtCol('ts', ['2026-01-01', '2026-01-04']), floatCol('v', [1, 2])];
    const entries = [
      entryFromDataFrame(makeDf('r1', cols())),
      entryFromDataFrame(makeDf('r2', cols())),
    ];
    const result = buildColumnComparison([targetFor(entries, 'ts')], entries,
      tsFor(entries, [undefined, undefined]))!;
    expect(result.timeSeriesUnit, 'days');
    expect(numbers(result.chartDf.getCol('ts (days)').toList()).join(','), '0,3,0,3');
  });

  test('partially configured targets fall back to the legacy axis', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [dtCol('ts', ['2026-01-01', '2026-01-02']), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('ts', [0, 60]), floatCol('v', [3, 4])])),
    ];
    const partial = new Map([[entries[0].id,
      new Map([[entries[0].nodes.tables[0].path, {mode: 'timeseries' as const}]])]]);
    const result = buildColumnComparison([targetFor(entries, 'ts')], entries, partial)!;
    expect(result.isKeyIndex, true);
    expect(result.timeSeriesUnit === undefined, true);
    expect(result.chartDf.getCol('ts').type, DG.COLUMN_TYPE.STRING);
  });

  test('null index cells are skipped for the start and stay null', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [null, 5, 10]), floatCol('v', [1, 2, 3])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [0, 5]), floatCol('v', [4, 5])])),
    ];
    const result = buildColumnComparison([targetFor(entries, 'time')], entries,
      tsFor(entries, ['s', 's']))!;
    const values = numbers(result.chartDf.getCol('time (s)').toList());
    expect(values[0], null);
    expect(values.slice(1).join(','), '0,5,0,5');
  });

  test('elapsed axis label avoids value column collisions', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [0, 1]), floatCol('time (s)', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [0, 1]), floatCol('time (s)', [3, 4])])),
    ];
    const result = buildColumnComparison([targetFor(entries, 'time')], entries,
      tsFor(entries, ['s', 's']))!;
    expect(result.indexColumnName, 'time (s) (2)');
    expect(numbers(result.chartDf.getCol('time (s) (2)').toList()).join(','), '0,1,0,1');
  });

  test('multi-value comparison shares the elapsed axis', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [10, 20]), floatCol('a', [1, 2]), floatCol('b', [5, 6])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [0, 10]), floatCol('a', [3, 4]), floatCol('b', [7, 8])])),
    ];
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, 'time'))
      .filter((t): t is ColumnTarget => t.kind === 'column');
    const result = buildColumnComparison(targets, entries, tsFor(entries, ['s', 's']))!;
    expect(result.timeSeriesUnit, 's');
    expect(numbers(result.chartDf.getCol('time (s)').toList()).join(','), '0,10,0,10');
    expect(result.valueColumnNames.join('|'), 'a|b');
  });
});

category('RunComparison: independent points dataframes', () => {
  const makeDf = (name: string, cols: DG.Column[]) => {
    const df = DG.DataFrame.fromColumns(cols);
    df.name = name;
    return df;
  };
  const floatCol = (name: string, values: (number | null)[]) =>
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
  const stringCol = (name: string, values: string[]) =>
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, values);
  const indexesFor = (entries: ComparisonEntry[], column: string) =>
    new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, column]])]));
  const splitsFor = (entries: ComparisonEntry[], column: string) =>
    new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, column]])]));
  const pointsFor = (entries: ComparisonEntry[]) =>
    new Map(entries.map((e) =>
      [e.id, new Map([[e.nodes.tables[0].path, {mode: 'points' as const}]])]));
  const columnTargetsFor = (
    entries: ComparisonEntry[], column: string, splits?: Map<string, Map<string, string>>,
  ) => matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries, column), splits)
    .filter((t): t is ColumnTarget => t.kind === 'column');

  test('all-points targets keep the raw index and set the flag', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [10, 20]), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [0, 10]), floatCol('v', [3, 4])])),
    ];
    const result = buildColumnComparison(
      [columnTargetsFor(entries, 'time')[0]], entries, pointsFor(entries))!;
    expect(result.isScatter, true);
    expect(result.isKeyIndex, false);
    expect(result.melted === undefined, true);
    expect(result.timeSeriesUnit === undefined, true);
    expect(result.indexColumnName, 'time');
    expect(result.chartDf.getCol('time').toList().join(','), '10,20,0,10');
    expect(result.chartDf.getCol('v').toList().join(','), '1,2,3,4');
  });

  test('partially configured targets chart as a line', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [10, 20]), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [0, 10]), floatCol('v', [3, 4])])),
    ];
    const partial = new Map([[entries[0].id,
      new Map([[entries[0].nodes.tables[0].path, {mode: 'points' as const}]])]]);
    const result = buildColumnComparison([columnTargetsFor(entries, 'time')[0]], entries, partial)!;
    expect(result.isScatter, false);
  });

  test('a key index disables points mode', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [stringCol('key', ['x', 'y']), floatCol('v', [1, 2])])),
      entryFromDataFrame(makeDf('r2', [stringCol('key', ['x', 'y']), floatCol('v', [3, 4])])),
    ];
    const result = buildColumnComparison(
      [columnTargetsFor(entries, 'key')[0]], entries, pointsFor(entries))!;
    expect(result.isKeyIndex, true);
    expect(result.isScatter, false);
  });

  test('multi-value points mode melts values into one column keyed by target name', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('time', [10, 20]), floatCol('a', [1, 2]), floatCol('b', [5, 6])])),
      entryFromDataFrame(makeDf('r2', [floatCol('time', [0, 10]), floatCol('a', [3, 4]), floatCol('b', [7, 8])])),
    ];
    const result = buildColumnComparison(
      columnTargetsFor(entries, 'time'), entries, pointsFor(entries))!;
    expect(result.isScatter, true);
    expect(result.melted!.seriesColumnName, 'Data');
    expect(result.melted!.valueColumnName, 'Value');
    expect(result.chartDf.rowCount, 8);
    expect(result.chartDf.getCol('time').toList().join(','), '10,20,0,10,10,20,0,10');
    expect(result.chartDf.getCol(RUN_COLUMN).toList().join(','), 'r1,r1,r2,r2,r1,r1,r2,r2');
    expect(result.chartDf.getCol('Data').toList().join(','), 'a,a,a,a,b,b,b,b');
    expect(result.chartDf.getCol('Value').toList().join(','), '1,2,3,4,5,6,7,8');
  });

  test('melted column names avoid collisions with the index', async () => {
    const entries = [
      entryFromDataFrame(makeDf('r1', [floatCol('Value', [0, 1]), floatCol('a', [1, 2]), floatCol('b', [5, 6])])),
      entryFromDataFrame(makeDf('r2', [floatCol('Value', [0, 1]), floatCol('a', [3, 4]), floatCol('b', [7, 8])])),
    ];
    const result = buildColumnComparison(
      columnTargetsFor(entries, 'Value'), entries, pointsFor(entries))!;
    expect(result.melted!.valueColumnName, 'Value (2)');
    expect(result.chartDf.getCol('Value (2)').toList().join(','), '1,2,3,4,5,6,7,8');
  });

  test('split columns carry into the melted frame', async () => {
    const cols = (values: number[]) => [
      floatCol('time', [0, 1]), stringCol('species', ['s1', 's2']),
      floatCol('a', values.slice(0, 2)), floatCol('b', values.slice(2)),
    ];
    const entries = [
      entryFromDataFrame(makeDf('r1', cols([1, 2, 5, 6]))),
      entryFromDataFrame(makeDf('r2', cols([3, 4, 7, 8]))),
    ];
    const result = buildColumnComparison(
      columnTargetsFor(entries, 'time', splitsFor(entries, 'species')), entries, pointsFor(entries))!;
    expect(result.splitColumnName, 'species');
    expect(result.chartDf.getCol('species').toList().join(','), 's1,s2,s1,s2,s1,s2,s1,s2');
  });
});

category('RunComparison: union of runs', () => {
  const col = (name: string, values: (number | null)[]) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
  const tableEntry = (name: string, columns: DG.Column[]) => {
    const df = DG.DataFrame.fromColumns(columns);
    df.name = name;
    return entryFromDataFrame(df);
  };
  const numbers = (values: any[]) => Array.from(values, (v) => v == null || isNaN(v) ? null : v);

  test('runs any selected column has are charted, whichever is selected first', async () => {
    const entries = [
      tableEntry('r1', [col('Time', [0, 1]), col('x', [1, 2]), col('y', [3, 4])]),
      tableEntry('r2', [col('Time', [0, 1]), col('x', [5, 6]), col('y', [7, 8])]),
      tableEntry('r3', [col('Time', [0, 1]), col('y', [9, 10])]),
    ];
    const indexes = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'Time']])]));
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    const x = targets.find((t) => t.displayName === 'x')!;
    const y = targets.find((t) => t.displayName === 'y')!;
    for (const order of [[x, y], [y, x]]) {
      const result = buildColumnComparison(order, entries)!;
      expectArray([...new Set(result.chartDf.getCol(RUN_COLUMN).toList())], ['r1', 'r2', 'r3']);
      expectArray(numbers(result.chartDf.getCol('x').toList()), [1, 2, 5, 6, null, null]);
    }
  });
});

category('RunComparison: unit conversion dataframes', () => {
  const col = (name: string, values: number[], units?: string) => {
    const column = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
    if (units)
      column.meta.units = units;
    return column;
  };
  const tableEntry = (name: string, columns: DG.Column[]) => {
    const df = DG.DataFrame.fromColumns(columns);
    df.name = name;
    return entryFromDataFrame(df);
  };
  const indexesFor = (entries: ComparisonEntry[]) =>
    new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'Time']])]));
  const scalar = (id: string, value: number, units: string) => {
    const e = scalarsEntry(id, id, []);
    e.nodes.scalars.push({path: 'p', name: 'p', valueType: 'double', units, value});
    e.nodes.scalars.push({path: 'q', name: 'q', valueType: 'double', value: 1});
    return e;
  };
  const convert = buildUnitConverter([
    {from: 'g/L', factor: 1000, to: 'mg/L'}, {from: 'bar', factor: 100, to: 'kPa'}, {from: 'h', factor: 60, to: 'min'},
  ]);

  test('a scalar converts into the first run\'s units', async () => {
    const entries = [scalar('a', 2, 'g/L'), scalar('b', 3000, 'mg/L')];
    const target = matchScalarTargets(entries.map((e) => e.nodes), undefined, convert)
      .find((t) => t.displayName === 'p')!;
    const result = buildScalarComparison(target, entries);
    expectArray(result.chartDf.getCol(result.valueColumnName).toList(), [2, 3]);
  });

  test('several scalars convert into their target units', async () => {
    const entries = [scalar('a', 1, 'bar'), scalar('b', 50, 'kPa')];
    const targets = matchScalarTargets(entries.map((e) => e.nodes), undefined, convert);
    const result = buildMultiScalarComparison(targets, entries)!;
    expectArray(result.chartDf.getCol('p (bar)').toList(), [1, 0.5]);
    expectArray(result.chartDf.getCol('q').toList(), [1, 1]);
  });

  test('column values convert into the first run\'s units, the index stays raw', async () => {
    const entries = [
      tableEntry('r1', [col('Time', [0, 1], 'h'), col('x', [0.5, 1], 'g/L')]),
      tableEntry('r2', [col('Time', [0, 30], 'min'), col('x', [250, 750], 'mg/L')]),
    ];
    const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries),
      undefined, undefined, undefined, convert);
    const result = buildColumnComparison([target], entries)!;
    expect(result.indexColumnName, 'Time');
    expect(result.valueColumnName, 'x (g/L)');
    expectArray(result.chartDf.getCol('Time').toList(), [0, 1, 0, 30]);
    expectArray(result.chartDf.getCol(result.valueColumnName).toList(), [0.5, 1, 0.25, 0.75]);
  });

  test('the melted scatter converts values', async () => {
    const entries = [
      tableEntry('r1', [col('Time', [0, 1]), col('x', [1, 2], 'bar'), col('y', [3, 4])]),
      tableEntry('r2', [col('Time', [0, 1]), col('x', [100, 200], 'kPa'), col('y', [5, 6])]),
    ];
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries),
      undefined, undefined, undefined, convert);
    const modes = new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, {mode: 'points' as const}]])]));
    const result = buildColumnComparison(targets, entries, modes)!;
    const melted = result.melted!;
    const values = result.chartDf.getCol(melted.valueColumnName).toList();
    const series = result.chartDf.getCol(melted.seriesColumnName).toList();
    expectArray(values.filter((_, i) => series[i] === 'x (bar)'), [1, 2, 1, 2]);
  });
});

category('RunComparison: units in column names', () => {
  const col = (name: string, values: number[], units?: string) => {
    const column = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
    if (units)
      column.meta.units = units;
    return column;
  };
  const tableEntry = (name: string, columns: DG.Column[]) => {
    const df = DG.DataFrame.fromColumns(columns);
    df.name = name;
    return entryFromDataFrame(df);
  };
  const indexesFor = (entries: ComparisonEntry[]) =>
    new Map(entries.map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'Time']])]));

  test('a scalar column carries its units, a unitless one stays plain', async () => {
    const entries = [scalarsEntry('a', 'a', []), scalarsEntry('b', 'b', [])];
    for (const e of entries)
      e.nodes.scalars.push({path: 'p', name: 'p', valueType: 'double', units: 'bar', value: 1});
    const [target] = matchScalarTargets(entries.map((e) => e.nodes));
    expect(buildScalarComparison(target, entries).valueColumnName, 'p (bar)');
    const plain = [scalarEntry('a', 'a', 1), scalarEntry('b', 'b', 2)];
    expect(buildScalarComparison(matchScalarTargets(plain.map((e) => e.nodes))[0], plain).valueColumnName, 'x');
  });

  test('value columns and a float index carry their units; a key index does not', async () => {
    const entries = [
      tableEntry('r1', [col('Time', [0, 1], 'h'), col('x', [1, 2], 'g/L'), col('y', [3, 4])]),
      tableEntry('r2', [col('Time', [0, 1], 'h'), col('x', [5, 6], 'g/L'), col('y', [7, 8])]),
    ];
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries));
    const result = buildColumnComparison(targets, entries)!;
    expect(result.indexColumnName, 'Time (h)');
    expectArray(result.valueColumnNames, ['x (g/L)', 'y']);
    expectArray(result.chartDf.columns.names(), ['Time (h)', RUN_COLUMN, 'x (g/L)', 'y']);
  });

  test('a float index carries units only when every run agrees on them', async () => {
    const build = (units: string[]) => {
      const entries = units.map((u, i) => tableEntry(`r${i}`, [col('Time', [0, 1], u), col('x', [1, 2])]));
      const [target] = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries));
      return buildColumnComparison([target], entries)!.indexColumnName;
    };
    expect(build(['h', 'hours']), 'Time (h)');
    expect(build(['h', 'min']), 'Time');
    expect(build(['h', '']), 'Time');
  });

  test('a units label never collides with another column', async () => {
    const entries = [
      tableEntry('r1', [col('Time', [0, 1], 'h'), col('x', [1, 2], 'g/L'), col('x (g/L)', [3, 4])]),
      tableEntry('r2', [col('Time', [0, 1], 'h'), col('x', [5, 6], 'g/L'), col('x (g/L)', [7, 8])]),
    ];
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexesFor(entries));
    const result = buildColumnComparison(targets, entries)!;
    expect(new Set(result.valueColumnNames).size, 2);
  });
});

category('RunComparison: comparison dispatch', () => {
  const scalarEntries = () => [
    scalarsEntry('a', 'r1', [{name: 'x', value: 1}, {name: 'y', value: 10}]),
    scalarsEntry('b', 'r2', [{name: 'x', value: 2}, {name: 'y', value: 20}]),
  ];
  const scalarTargets = (entries: ComparisonEntry[]) => {
    const targets = matchScalarTargets(entries.map((e) => e.nodes));
    return {x: targets.find((t) => t.displayName === 'x')!, y: targets.find((t) => t.displayName === 'y')!};
  };
  const columnEntries = () => [1, 2].map((run) => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Time', [0, 1]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [run, run + 1]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [run * 10, run * 10 + 1]),
    ]);
    df.name = `r${run}`;
    return entryFromDataFrame(df);
  });
  const columnTargets = (entries: ComparisonEntry[]) => {
    const indexes = new Map(entries.filter((e) => e.nodes.tables.length > 0)
      .map((e) => [e.id, new Map([[e.nodes.tables[0].path, 'Time']])]));
    const targets = matchColumnTargets(entries.map((e) => e.nodes), indexes);
    return {x: targets.find((t) => t.displayName === 'x')!, y: targets.find((t) => t.displayName === 'y')!};
  };

  test('no selection or a single run builds nothing', async () => {
    const entries = scalarEntries();
    expect(buildComparison([], entries), null);
    expect(buildComparison([scalarTargets(entries).x], entries.slice(0, 1)), null);
  });

  test('one scalar builds a scalar comparison', async () => {
    const entries = scalarEntries();
    const result = buildComparison([scalarTargets(entries).x], entries)!;
    expect(result.kind, 'scalar');
    expect(result.chartDf.name, 'Comparison: x');
    expect(result.chartDf.col('Path') != null, true);
  });

  test('several scalars build a multi-scalar comparison, anchor first', async () => {
    const entries = scalarEntries();
    const {x, y} = scalarTargets(entries);
    const result = buildComparison([y, x], entries)!;
    if (result.kind !== 'multi-scalar')
      throw new Error(`expected multi-scalar, got ${result.kind}`);
    expectArray(result.valueColumnNames, ['y', 'x']);
  });

  test('one column builds a single-value column comparison', async () => {
    const entries = columnEntries();
    const result = buildComparison([columnTargets(entries).x], entries)!;
    if (result.kind !== 'column')
      throw new Error(`expected column, got ${result.kind}`);
    expectArray(result.valueColumnNames, ['x']);
    expect(result.valueColumnName, 'x');
    expect(result.chartDf.name, 'Comparison: x');
  });

  test('several columns build one multi-value column comparison', async () => {
    const entries = columnEntries();
    const {x, y} = columnTargets(entries);
    const result = buildComparison([y, x], entries)!;
    if (result.kind !== 'column')
      throw new Error(`expected column, got ${result.kind}`);
    expectArray(result.valueColumnNames, ['y', 'x']);
    expect(result.chartDf.name, 'Comparison: multiple values');
  });

  test('targets of another kind than the anchor are ignored', async () => {
    const entries = [...scalarEntries(), ...columnEntries()];
    const scalar = scalarTargets(entries).x;
    const column = columnTargets(entries).x;
    expect(buildComparison([scalar, column], entries)!.kind, 'scalar');
    const result = buildComparison([column, scalar], entries)!;
    if (result.kind !== 'column')
      throw new Error(`expected column, got ${result.kind}`);
    expectArray(result.valueColumnNames, ['x']);
  });
});

category('RunComparison: column primitives', () => {
  test('missing cells are null, datetimes epoch ms, bigints and strings as is', async () => {
    const day = dayjs('2026-01-02T03:04:05Z');
    const floats = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'f', [1.5, null]);
    const dates = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'd', [day, null]);
    const bigints = DG.Column.fromBigInt64Array('b', new BigInt64Array([BigInt(7)]));
    const strings = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 's', ['x', '']);
    const [f, d, b, s] = [floats, dates, bigints, strings].map(columnPrimitives);
    expect(f[0], 1.5);
    expect(f[1] === null, true);
    expect(d[0], day.valueOf());
    expect(d[1] === null, true);
    expect(b[0] === BigInt(7), true);
    expect(s[0], 'x');
    expect(s[1] === null, true);
  });
});
