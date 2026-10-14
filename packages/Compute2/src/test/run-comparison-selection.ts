import {category, test, expect, expectArray, expectObject} from '@datagrok-libraries/test/src/test';
import {
  matchScalarTargets, matchColumnTargets, buildUnitConverter, buildAliasGroups,
} from '../components/RunComparison/matching';
import {
  getEntryStatuses, matchesFilter, compatibleTargetsFor,
  isSplitCandidate, selectionToMap, computeIndexRows, isTargetEqualAcrossRuns,
  parseTimeUnit, resolveAxisModes, targetAxisMode, resolveDisplayUnit, pickAutoUnit, targetContexts,
  selectionGaps, withAnnotatedDefaults, PickerSelections, unmatchedItems,
} from '../components/RunComparison/selection';
import {ColumnBinding, ColumnTarget, AxisConfigMap, TimeUnit, TableNodeInfo} from '../components/RunComparison/types';
import {makeEntry, indexMap} from './run-comparison-fixtures';

category('RunComparison: entry statuses', () => {
  test('reports matched, index-not-set, and no-similar-data', async () => {
    const withTable = {path: 't', columns: [{name: 'time'}, {name: 'value'}]};
    const entries = [
      makeEntry('a', [], [withTable]),
      makeEntry('b', [], [withTable]),
      makeEntry('c', [], [{path: 'other', columns: [{name: 'x'}, {name: 'y'}]}]),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [target] = matchColumnTargets(entries, indexes);
    const statuses = getEntryStatuses(entries, target, indexes);
    expect(statuses.find((s) => s.entryId === 'a')!.matched, true);
    expect(statuses.find((s) => s.entryId === 'b')!.matched, true);
    const statusC = statuses.find((s) => s.entryId === 'c')!;
    expect(statusC.matched, false);
    expect(statusC.reason, 'index not set');
  });

  test('reports no similar data when target is absent', async () => {
    const entries = [makeEntry('a', [{name: 'x'}])];
    const statuses = getEntryStatuses(entries, null, indexMap({}));
    expect(statuses[0].matched, false);
    expect(statuses[0].reason, 'no similar data');
  });

  test('index set but no similar column reports no-similar-data', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]}]),
      makeEntry('b', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]}]),
      makeEntry('c', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'unrelated'}]}]),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}, c: {t: 'time'}});
    const [target] = matchColumnTargets(entries, indexes);
    const statusC = getEntryStatuses(entries, target, indexes).find((s) => s.entryId === 'c')!;
    expect(statusC.matched, false);
    expect(statusC.reason, 'no similar data');
  });

  test('entry with all candidates toggled off reports disabled', async () => {
    const table = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [base] = matchColumnTargets(entries, indexes);
    const [target] = matchColumnTargets(entries, indexes, undefined,
      {[base.key]: {'b|t|height': false}});
    const statusB = getEntryStatuses(entries, target, indexes).find((s) => s.entryId === 'b')!;
    expect(statusB.matched, false);
    expect(statusB.reason, 'disabled');
  });

  test('a column target on an entry without tables reports no similar data', async () => {
    const table = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table]), makeEntry('c', [{name: 'x'}])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [target] = matchColumnTargets(entries, indexes);
    expect(getEntryStatuses(entries, target, indexes).find((s) => s.entryId === 'c')!.reason, 'no similar data');
  });

  test('a disabled run reports disabled even when another table has the value in other units', async () => {
    const table = (units: string) => ({path: units, columns: [{name: 'time', type: 'int'}, {name: 'height', units}]});
    const entries = [makeEntry('a', [], [table('m')]), makeEntry('b', [], [table('m'), table('ft')])];
    const indexes = indexMap({a: {m: 'time'}, b: {m: 'time', ft: 'time'}});
    const [base] = matchColumnTargets(entries, indexes);
    const [target] = matchColumnTargets(entries, indexes, undefined, {[base.key]: {'b|m|height': false}});
    expect(getEntryStatuses(entries, target, indexes).find((s) => s.entryId === 'b')!.reason, 'disabled');
  });

  test('scalar target unmatched entry reports no-similar-data', async () => {
    const entries = [
      makeEntry('a', [{name: 'x'}]),
      makeEntry('b', [{name: 'x'}]),
      makeEntry('c', [{name: 'other'}]),
    ];
    const [target] = matchScalarTargets(entries);
    const statusC = getEntryStatuses(entries, target, indexMap({})).find((s) => s.entryId === 'c')!;
    expect(statusC.matched, false);
    expect(statusC.reason, 'no similar data');
  });
});

category('RunComparison: filters and compatibility', () => {
  test('matchesFilter: empty, substring, fuzzy token, miss', async () => {
    expect(matchesFilter('', 'anything'), true);
    expect(matchesFilter('temp', 'Init_Temperature'), true);
    expect(matchesFilter('temperatures', 'Init Temperature'), true);
    expect(matchesFilter('pressure', 'height'), false);
  });

  const bind = (entryId: string, tablePath: string): ColumnBinding =>
    ({entryId, tablePath, tableName: tablePath, columnName: 'v', indexColumnName: 'time'});
  const columnTarget = (key: string, bindings: ColumnBinding[]): ColumnTarget => ({
    kind: 'column', key, displayName: key, confidence: 'exact', unitsWarning: false,
    coverage: bindings.length, total: 3,
    candidates: [], bindings,
  });

  test('columns without shared runs are suggested like scalars', async () => {
    const anchor = columnTarget('anchor', [bind('a', 't')]);
    const other = columnTarget('other', [bind('b', 't')]);
    expect(compatibleTargetsFor(anchor, [anchor, other], () => 'int').length, 2);
  });

  test('a column with a non-line-chartable index is not suggested', async () => {
    const anchor = columnTarget('anchor', [bind('a', 't'), bind('b', 't')]);
    const other = columnTarget('other', [bind('a', 't'), bind('b', 'keyed')]);
    const getType = (_entryId: string, tablePath: string) => tablePath === 'keyed' ? 'string' : 'int';
    expectArray(compatibleTargetsFor(anchor, [anchor, other], getType).map((t) => t.key), ['anchor']);
  });

  test('cross-table pick stays suggested', async () => {
    const full = {path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'height'}, {name: 'velocity'},
    ]};
    const extra = {path: 'u', columns: [{name: 'time', type: 'int'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [full]), makeEntry('b', [], [full, extra])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time', u: 'time'}});
    const base = matchColumnTargets(entries, indexes);
    const heightKey = base.find((t) => t.displayName === 'height')!.key;
    const edited = matchColumnTargets(entries, indexes, undefined,
      {[heightKey]: {'b|u|height': true}});
    const anchor = edited.find((t) => t.displayName === 'height')!;
    const velocity = edited.find((t) => t.displayName === 'velocity')!;
    expect(anchor.bindings.find((b) => b.entryId === 'b')!.tablePath, 'u');
    expect(velocity.bindings.find((b) => b.entryId === 'b')!.tablePath, 't');
    expect(compatibleTargetsFor(anchor, edited, () => 'int').length, 2);
  });

  test('a column anchor suggests every line-indexed column', async () => {
    const table = {path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'height'}, {name: 'velocity'}, {name: 'other'},
    ]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const targets = matchColumnTargets(entries, indexes);
    const anchor = targets.find((t) => t.displayName === 'height')!;
    const getType = () => 'int';
    const compatible = compatibleTargetsFor(anchor, targets, getType);
    expect(compatible.length, 3);
    expect(compatibleTargetsFor(null, targets, getType).length, 0);
  });

  test('compatibleTargetsFor rejects non-line-chartable index', async () => {
    const table = {path: 't', columns: [{name: 'key', type: 'string'}, {name: 'height'}, {name: 'velocity'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const indexes = indexMap({a: {t: 'key'}, b: {t: 'key'}});
    const targets = matchColumnTargets(entries, indexes);
    expect(compatibleTargetsFor(targets[0], targets, () => 'string').length, 0);
    expect(compatibleTargetsFor(targets[0], targets, () => 'datetime').length, 2);
  });

  test('raw fuzzy siblings do not block multi-value mode', async () => {
    const wfTable = {path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'temperature'}, {name: 'velocity'},
    ]};
    const rawTable = {path: 'exp', columns: [
      {name: 'time', type: 'int'}, {name: 'temperature'}, {name: 'temperatures'},
    ]};
    const entries = [
      makeEntry('a', [], [wfTable]),
      makeEntry('b', [], [wfTable]),
      makeEntry('raw', [], [rawTable], 'raw'),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}, raw: {exp: 'time'}});
    const targets = matchColumnTargets(entries, indexes);
    const anchor = targets.find((t) => t.displayName === 'temperature')!;
    expect(anchor.bindings.filter((b) => b.entryId === 'raw').length, 1);
    expect(compatibleTargetsFor(anchor, targets, () => 'int').length >= 1, true);
  });

  test('scalar anchor suggests all scalars and no columns', async () => {
    const table = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]};
    const entries = [
      makeEntry('a', [{name: 'x'}, {name: 'y'}], [table]),
      makeEntry('b', [{name: 'x'}, {name: 'y'}], [table]),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const targets = [...matchScalarTargets(entries), ...matchColumnTargets(entries, indexes)];
    const anchor = targets.find((t) => t.kind === 'scalar')!;
    const compatible = compatibleTargetsFor(anchor, targets, () => 'int');
    expect(compatible.length, 2);
    expect(compatible.every((t) => t.kind === 'scalar'), true);
  });

  test('column anchor never suggests scalars', async () => {
    const table = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}, {name: 'velocity'}]};
    const entries = [
      makeEntry('a', [{name: 'x'}], [table]),
      makeEntry('b', [{name: 'x'}], [table]),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const targets = [...matchScalarTargets(entries), ...matchColumnTargets(entries, indexes)];
    const anchor = targets.find((t) => t.kind === 'column')!;
    const compatible = compatibleTargetsFor(anchor, targets, () => 'int');
    expect(compatible.length, 2);
    expect(compatible.every((t) => t.kind === 'column'), true);
  });

  test('isSplitCandidate requires a string column that is not the index', async () => {
    expect(isSplitCandidate({name: 'species', type: 'string'}, 'time'), true);
    expect(isSplitCandidate({name: 'time', type: 'string'}, 'time'), false);
    expect(isSplitCandidate({name: 'height', type: 'double'}, 'time'), false);
  });
});

category('RunComparison: index rows', () => {
  const columns = [{name: 'time', type: 'int'}, {name: 'height'}, {name: 'species', type: 'string'}];

  test('merging groups same-function outputs and intersects columns', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim',
        columns: [...columns, {name: 'extra', type: 'int'}]}]),
    ];
    const rows = computeIndexRows(entries, {}, {}, true, {});
    expect(rows.length, 1);
    expect(rows[0].coverage!.count, 2);
    expect(rows[0].members.length, 2);
    expect(rows[0].candidates.some((c) => c.name === 'extra'), false);
  });

  test('merged rows of tables at different paths use the table name and count runs', async () => {
    const tables = [
      {path: 'days/day[1]/df', name: 'Profile', nqName: 'Pkg:Day', columns},
      {path: 'days/day[2]/df', name: 'Profile', nqName: 'Pkg:Day', columns},
    ];
    const entries = [makeEntry('a', [], tables), makeEntry('b', [], tables.slice(0, 1))];
    const rows = computeIndexRows(entries, {}, {}, true, {});
    expect(rows.length, 1);
    expect(rows[0].label, 'Profile');
    expect(rows[0].members.length, 3);
    expect(rows[0].coverage!.count, 2);
    expect(rows[0].coverage!.total, 2);
  });

  test('tables without nqName never merge', async () => {
    const entries = [
      makeEntry('a', [], [{path: 'raw', columns}]),
      makeEntry('b', [], [{path: 'raw', columns}]),
    ];
    const rows = computeIndexRows(entries, {}, {}, true, {});
    expect(rows.length, 2);
    expect(rows[0].coverage === undefined, true);
  });

  test('merge off yields one row per table', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
    ];
    expect(computeIndexRows(entries, {}, {}, false, {}).length, 2);
  });

  test('stale and mixed selections fall back to unset', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
    ];
    const stale = computeIndexRows(entries, {a: {'s/df': 'removed'}, b: {'s/df': 'removed'}}, {}, true, {});
    expect(stale[0].current, '');
    const mixed = computeIndexRows(entries, {a: {'s/df': 'time'}, b: {'s/df': 'height'}}, {}, true, {});
    expect(mixed[0].current, '');
    const agreed = computeIndexRows(entries, {a: {'s/df': 'time'}, b: {'s/df': 'time'}}, {}, true, {});
    expect(agreed[0].current, 'time');
  });

  test('split candidates exclude the chosen index and non-strings', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
    ];
    const rows = computeIndexRows(entries, {a: {'s/df': 'time'}, b: {'s/df': 'time'}},
      {a: {'s/df': 'species'}, b: {'s/df': 'species'}}, true, {});
    expectArray(rows[0].splitCandidates.map((c) => c.name), ['species']);
    expect(rows[0].currentSplit, 'species');
  });

  test('numeric, datetime, and string columns are always offered as index candidates', async () => {
    const table = {path: 's/df', nqName: 'Pkg:Sim', columns: [
      {name: 'time_f'}, {name: 'step', type: 'int'}, {name: 'stamp', type: 'datetime'},
      {name: 'label', type: 'string'}, {name: 'flag', type: 'bool'},
    ]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const rows = computeIndexRows(entries, {}, {}, true, {});
    expectArray(rows[0].candidates.map((c) => c.name), ['time_f', 'step', 'stamp', 'label']);
  });

  test('annotated default index is offered even for an exotic type', async () => {
    const table = {path: 's/df', nqName: 'Pkg:Sim', defaultIndexColumn: 'flag',
      columns: [{name: 'flag', type: 'bool'}, {name: 'step', type: 'int'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const rows = computeIndexRows(entries, {}, {}, true, {});
    expectArray(rows[0].candidates.map((c) => c.name), ['flag', 'step', 'height']);
  });

  test('selectionToMap keeps only valid selections', async () => {
    const map = selectionToMap({a: {t: 'time', u: 'bad'}, b: {t: ''}}, (_e, _t, col) => col !== 'bad');
    expect(map.get('a')!.get('t'), 'time');
    expect(map.get('a')!.has('u'), false);
    expect(map.get('b')!.size, 0);
  });
});

category('RunComparison: equal values', () => {
  const scalarTarget = (values: (number | null)[]) => {
    const entries = values.map((value, i) => makeEntry(`e${i}`, [{name: 'x', value}]));
    return matchScalarTargets(entries)[0];
  };

  const binding = (
    entryId: string, columnName = 'height', splitColumnName?: string,
  ): ColumnBinding => ({
    entryId, tablePath: 't', tableName: 't', columnName, indexColumnName: 'time', splitColumnName,
  });
  const columnTarget = (bindings: ColumnBinding[]): ColumnTarget => ({
    kind: 'column', key: 'k', displayName: 'height', confidence: 'exact', unitsWarning: false,
    coverage: bindings.length, total: bindings.length,
    candidates: [], bindings,
  });
  const getter = (data: Record<string, Record<string, unknown[]>>) =>
    (entryId: string, _tablePath: string, columnName: string) => data[entryId]?.[columnName];

  test('scalar target with identical values is equal', async () => {
    expect(isTargetEqualAcrossRuns(scalarTarget([5, 5, 5]), () => undefined), true);
    expect(isTargetEqualAcrossRuns(scalarTarget([5, 6]), () => undefined), false);
    expect(isTargetEqualAcrossRuns(scalarTarget([null, null]), () => undefined), true);
  });

  test('numbers within 0.1% relative tolerance are equal', async () => {
    expect(isTargetEqualAcrossRuns(scalarTarget([1000, 1000.5]), () => undefined), true);
    expect(isTargetEqualAcrossRuns(scalarTarget([1000, 1002]), () => undefined), false);
    expect(isTargetEqualAcrossRuns(scalarTarget([-1000, -1000.5]), () => undefined), true);
    expect(isTargetEqualAcrossRuns(scalarTarget([0, 1e-9]), () => undefined), false);
    expect(isTargetEqualAcrossRuns(scalarTarget([5, null]), () => undefined), false);
    const target = columnTarget([binding('a'), binding('b')]);
    expect(isTargetEqualAcrossRuns(target, getter({
      a: {height: [1000, 2000], time: [0, 1]}, b: {height: [1000.5, 1999], time: [0, 1]},
    })), true);
    expect(isTargetEqualAcrossRuns(target, getter({
      a: {height: [1000, 2000], time: [0, 1]}, b: {height: [1000.5, 1990], time: [0, 1]},
    })), false);
  });

  test('single-run target is never equal', async () => {
    const target = columnTarget([binding('a')]);
    expect(isTargetEqualAcrossRuns(target, getter({a: {height: [1], time: [0]}})), false);
  });

  test('column target with identical value and index data is equal', async () => {
    const target = columnTarget([binding('a'), binding('b', 'height2')]);
    const data = {
      a: {height: [1, 2, null], time: [0, 1, 2]},
      b: {height2: [1, 2, null], time: [0, 1, 2]},
    };
    expect(isTargetEqualAcrossRuns(target, getter(data)), true);
  });

  test('differing values, index, or length break equality', async () => {
    const target = columnTarget([binding('a'), binding('b')]);
    expect(isTargetEqualAcrossRuns(target, getter({
      a: {height: [1, 2], time: [0, 1]}, b: {height: [1, 3], time: [0, 1]},
    })), false);
    expect(isTargetEqualAcrossRuns(target, getter({
      a: {height: [1, 2], time: [0, 1]}, b: {height: [1, 2], time: [0, 2]},
    })), false);
    expect(isTargetEqualAcrossRuns(target, getter({
      a: {height: [1, 2], time: [0, 1]}, b: {height: [1], time: [0]},
    })), false);
  });

  test('unfetchable column data counts as differing', async () => {
    const target = columnTarget([binding('a'), binding('b')]);
    expect(isTargetEqualAcrossRuns(target, getter({a: {height: [1], time: [0]}})), false);
  });

  test('datetime indexes compare exactly, not relatively', async () => {
    const target = columnTarget([binding('a'), binding('b')]);
    const day = 86400000;
    const t0 = Date.UTC(2026, 7, 1);
    const shifted = {
      a: {height: [1, 2], time: [t0, t0 + 60000]},
      b: {height: [1, 2], time: [t0 + day, t0 + day + 60000]},
    };
    const typed = (_entryId: string, _tablePath: string, columnName: string) =>
      columnName === 'time' ? 'datetime' : 'double';
    // a one-day shift is within 0.1% of an epoch-ms timestamp — exact compare must catch it
    expect(isTargetEqualAcrossRuns(target, getter(shifted), typed), false);
    expect(isTargetEqualAcrossRuns(target, getter({a: shifted.a, b: shifted.a}), typed), true);
    // numeric indexes keep the relative tolerance
    expect(isTargetEqualAcrossRuns(target, getter({
      a: {height: [1, 2], time: [1e6, 1e6 + 1]}, b: {height: [1, 2], time: [1e6 + 100, 1e6 + 101]},
    }), () => 'double'), true);
  });

  test('split columns must match on both sides', async () => {
    const withSplit = columnTarget([binding('a', 'height', 'species'), binding('b', 'height', 'species')]);
    const same = {height: [1, 2], time: [0, 1], species: ['s1', 's2']};
    expect(isTargetEqualAcrossRuns(withSplit, getter({a: same, b: same})), true);
    expect(isTargetEqualAcrossRuns(withSplit, getter({
      a: same, b: {...same, species: ['s1', 's3']},
    })), false);
    const mixed = columnTarget([binding('a', 'height', 'species'), binding('b')]);
    expect(isTargetEqualAcrossRuns(mixed, getter({a: same, b: same})), false);
  });
});

category('RunComparison: axis mode selection', () => {
  const floatTable = (units?: string) => ({path: 't', columns: [
    {name: 'time', type: 'double', units}, {name: 'height'},
  ]});
  const dtTable = {path: 't', columns: [{name: 'stamp', type: 'datetime'}, {name: 'height'}]};

  const bind = (entryId: string, tablePath = 't'): ColumnBinding =>
    ({entryId, tablePath, tableName: tablePath, columnName: 'height', indexColumnName: 'time'});
  const target = (bindings: ColumnBinding[]): ColumnTarget => ({
    kind: 'column', key: 'k', displayName: 'height', confidence: 'exact', unitsWarning: false,
    coverage: bindings.length, total: bindings.length,
    candidates: [], bindings,
  });
  const config = (
    spec: Record<string, Record<string, {mode?: 'timeseries' | 'points', units?: TimeUnit}>>,
  ): AxisConfigMap =>
    new Map(Object.entries(spec).map(([entryId, tables]) => [entryId,
      new Map(Object.entries(tables).map(([path, {mode, units}]) =>
        [path, {mode: mode ?? 'timeseries', units}]))]));

  test('parseTimeUnit recognizes aliases and rejects ambiguity', async () => {
    expect(parseTimeUnit('s'), 's');
    expect(parseTimeUnit(' SEC '), 's');
    expect(parseTimeUnit('minutes'), 'min');
    expect(parseTimeUnit('hr'), 'h');
    expect(parseTimeUnit('day'), 'days');
    expect(parseTimeUnit('ms'), 'ms');
    expect(parseTimeUnit('Hours'), 'h');
    expect(parseTimeUnit('m') === undefined, true);
    expect(parseTimeUnit('H') === undefined, true);
    expect(parseTimeUnit('meters') === undefined, true);
    expect(parseTimeUnit(undefined) === undefined, true);
  });

  test('index rows carry axis controls for numeric and datetime indexes', async () => {
    const entries = [makeEntry('a', [], [floatTable()]), makeEntry('b', [], [dtTable])];
    const rows = computeIndexRows(entries, {a: {t: 'time'}, b: {t: 'stamp'}}, {}, false, {});
    const floatRow = rows.find((row) => row.entryName === 'a')!;
    expect(floatRow.axis!.mode, 'series');
    expect(floatRow.axis!.units, 's');
    const dtRow = rows.find((row) => row.entryName === 'b')!;
    expect(dtRow.axis!.mode, 'series');
    expect(dtRow.axis!.units === undefined, true);
  });

  test('no controls for string or unset indexes', async () => {
    const table = {path: 't', columns: [{name: 'key', type: 'string'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [floatTable()])];
    const rows = computeIndexRows(entries, {a: {t: 'key'}}, {}, false, {});
    expect(rows.find((row) => row.entryName === 'a')!.axis === undefined, true);
    expect(rows.find((row) => row.entryName === 'b')!.axis === undefined, true);
  });

  test('units prefill from column metadata, stored pick wins', async () => {
    const entries = [makeEntry('a', [], [floatTable('hours')]), makeEntry('b', [], [floatTable('sec')])];
    const selection = {a: {t: 'time'}, b: {t: 'time'}};
    const rows = computeIndexRows(entries, selection, {}, false, {});
    expect(rows.find((row) => row.entryName === 'a')!.axis!.units, 'h');
    expect(rows.find((row) => row.entryName === 'b')!.axis!.units, 's');
    const stored = computeIndexRows(entries, selection, {}, false,
      {a: {t: {mode: 'timeseries', units: 'min'}}});
    const rowA = stored.find((row) => row.entryName === 'a')!;
    expect(rowA.axis!.mode, 'timeseries');
    expect(rowA.axis!.units, 'min');
  });

  test('stale config is ignored on a non-time index and resurfaces later', async () => {
    const table = {path: 't', columns: [
      {name: 'time'}, {name: 'key', type: 'string'}, {name: 'height'},
    ]};
    const entries = [makeEntry('a', [], [table])];
    const stored = {a: {t: {mode: 'points' as const, units: 'min' as const}}};
    const onString = computeIndexRows(entries, {a: {t: 'key'}}, {}, false, stored);
    expect(onString[0].axis === undefined, true);
    const back = computeIndexRows(entries, {a: {t: 'time'}}, {}, false, stored);
    expect(back[0].axis!.mode, 'points');
    expect(back[0].axis!.units, 'min');
  });

  test('elapsed on a datetime index stores no units, so the prefill wins after an index switch', async () => {
    const table = {path: 't', columns: [
      {name: 'stamp', type: 'datetime'}, {name: 'time', type: 'double', units: 'hours'}, {name: 'height'},
    ]};
    const entries = [makeEntry('a', [], [table])];
    const stored = {a: {t: {mode: 'timeseries' as const}}};
    const onDatetime = computeIndexRows(entries, {a: {t: 'stamp'}}, {}, false, stored);
    expect(onDatetime[0].axis!.mode, 'timeseries');
    expect(onDatetime[0].axis!.units === undefined, true);
    const onNumeric = computeIndexRows(entries, {a: {t: 'time'}}, {}, false, stored);
    expect(onNumeric[0].axis!.units, 'h');
    const resolved = resolveAxisModes(entries, indexMap({a: {t: 'time'}}), stored);
    expect(resolved.get('a')!.get('t')!.units, 'h');
  });

  test('merged rows agree on the mode like index picks', async () => {
    const table = {path: 's/df', nqName: 'Pkg:Sim', columns: [{name: 'time'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const selection = {a: {'s/df': 'time'}, b: {'s/df': 'time'}};
    const both = computeIndexRows(entries, selection, {}, true, {
      a: {'s/df': {mode: 'points'}}, b: {'s/df': {mode: 'points'}},
    });
    expect(both[0].axis!.mode, 'points');
    const one = computeIndexRows(entries, selection, {}, true, {
      a: {'s/df': {mode: 'points'}},
    });
    expect(one[0].axis!.mode, 'series');
    const mixed = computeIndexRows(entries, selection, {}, true, {
      a: {'s/df': {mode: 'points'}}, b: {'s/df': {mode: 'timeseries', units: 's'}},
    });
    expect(mixed[0].axis!.mode, 'series');
  });

  test('resolveAxisModes keeps only non-series tables with a time-typed index', async () => {
    const entries = [
      makeEntry('a', [], [floatTable()]),
      makeEntry('b', [], [dtTable]),
      makeEntry('c', [], [floatTable()]),
      makeEntry('d', [], [floatTable()]),
      makeEntry('e', [], [floatTable('sec')]),
    ];
    const map = resolveAxisModes(entries,
      indexMap({a: {t: 'time'}, b: {t: 'stamp'}, d: {t: 'time'}, e: {t: 'time'}}), {
        a: {t: {mode: 'timeseries', units: 'min'}},
        b: {t: {mode: 'timeseries', units: 's'}},
        c: {t: {mode: 'timeseries', units: 's'}},
        d: {t: {mode: 'series', units: 's'}},
        e: {t: {mode: 'points'}},
      });
    expect(map.get('a')!.get('t')!.mode, 'timeseries');
    expect(map.get('a')!.get('t')!.units, 'min');
    expect(map.get('b')!.get('t')!.units === undefined, true);
    expect(map.has('c'), false);
    expect(map.has('d'), false);
    // scatter resolves without units even on a numeric index with metadata
    expect(map.get('e')!.get('t')!.mode, 'points');
    expect(map.get('e')!.get('t')!.units === undefined, true);
  });

  test('targetAxisMode reports full, partial, none per mode', async () => {
    const tgt = target([bind('a'), bind('b')]);
    expect(targetAxisMode(tgt, config({a: {t: {units: 's'}}, b: {t: {}}}), 'timeseries'), 'full');
    expect(targetAxisMode(tgt, config({a: {t: {units: 's'}}}), 'timeseries'), 'partial');
    expect(targetAxisMode(tgt, config({}), 'timeseries'), 'none');
    expect(targetAxisMode(target([]), config({a: {t: {units: 's'}}}), 'timeseries'), 'none');
    const scatter = config({a: {t: {mode: 'points'}}, b: {t: {mode: 'points'}}});
    expect(targetAxisMode(tgt, scatter, 'points'), 'full');
    expect(targetAxisMode(tgt, scatter, 'timeseries'), 'none');
    const mixed = config({a: {t: {mode: 'points'}}, b: {t: {units: 's'}}});
    expect(targetAxisMode(tgt, mixed, 'points'), 'partial');
    expect(targetAxisMode(tgt, mixed, 'timeseries'), 'partial');
  });

  test('resolveDisplayUnit prefers the first numeric elapsed binding, datetime-only is auto', async () => {
    const bindings = [bind('a'), bind('b')];
    expect(resolveDisplayUnit(bindings, config({a: {t: {}}, b: {t: {units: 'min'}}})), 'min');
    expect(resolveDisplayUnit(bindings, config({a: {t: {units: 'h'}}, b: {t: {units: 'min'}}})), 'h');
    expect(resolveDisplayUnit(bindings, config({a: {t: {}}, b: {t: {}}})), 'auto');
    expect(resolveDisplayUnit(bindings,
      config({a: {t: {mode: 'points', units: 'min'}}, b: {t: {units: 'h'}}})), 'h');
  });

  test('pickAutoUnit boundaries', async () => {
    expect(pickAutoUnit(3 * 86_400_000), 'days');
    expect(pickAutoUnit(2 * 60_000), 'min');
    expect(pickAutoUnit(2 * 60_000 - 1), 's');
    expect(pickAutoUnit(2000), 's');
    expect(pickAutoUnit(500), 'ms');
  });

  test('partial targets flag unconfigured matched runs with a warning', async () => {
    const withTable = {path: 't', columns: [{name: 'time'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [withTable]), makeEntry('b', [], [withTable])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [tgt] = matchColumnTargets(entries, indexes);
    const partial = getEntryStatuses(entries, tgt, indexes, {axisModes: config({a: {t: {units: 's'}}})});
    expect(partial.find((s) => s.entryId === 'a')!.warning === undefined, true);
    expect(partial.find((s) => s.entryId === 'b')!.warning, 'relative timeseries not set');
    const partialScatter = getEntryStatuses(entries, tgt, indexes, {axisModes: config({a: {t: {mode: 'points'}}})});
    expect(partialScatter.find((s) => s.entryId === 'b')!.warning, 'independent points not set');
    const full = getEntryStatuses(entries, tgt, indexes,
      {axisModes: config({a: {t: {units: 's'}}, b: {t: {units: 's'}}})});
    expect(full.every((s) => s.warning === undefined), true);
    const none = getEntryStatuses(entries, tgt, indexes, {axisModes: config({})});
    expect(none.every((s) => s.warning === undefined), true);
  });
});

category('RunComparison: step contexts', () => {
  test('repeated names get their step as context, unique names none', async () => {
    const day = (n: number) =>
      ({name: 'Final biomass', friendlyPath: `Daily cultivation · Day ${n} · Final biomass`});
    const targets = matchScalarTargets([
      makeEntry('a', [day(1), day(2), {name: 'Volume'}]), makeEntry('b', [day(1), day(2), {name: 'Volume'}]),
    ]);
    const contexts = targetContexts(targets);
    expectArray(targets.filter((t) => t.displayName === 'Final biomass').map((t) => contexts.get(t.key)!),
      ['Daily cultivation · Day 1', 'Daily cultivation · Day 2']);
    expect(contexts.has(targets.find((t) => t.displayName === 'Volume')!.key), false);
  });
});

category('RunComparison: unmatched items', () => {
  const unmatched = (entries: ReturnType<typeof makeEntry>[], indexes = indexMap({})) => unmatchedItems(entries,
    [...matchScalarTargets(entries), ...matchColumnTargets(entries, indexes)], indexes);

  test('a scalar in only one run is unmatched', async () => {
    const items = unmatched([makeEntry('a', [{name: 'x'}, {name: 'y'}]), makeEntry('b', [{name: 'x'}])]);
    expect(items.length, 1);
    expectObject(items[0], {name: 'y', runCount: 1});
    expect(items[0].units === undefined, true);
  });

  test('a column in other units is unmatched with its units', async () => {
    const table = (units: string) => ({path: 't', columns: [{name: 'time'}, {name: 'Biomass', units}]});
    const items = unmatched(
      [makeEntry('a', [], [table('g/L')]), makeEntry('b', [], [table('g/L')]), makeEntry('c', [], [table('mg/L')])],
      indexMap({a: {t: 'time'}, b: {t: 'time'}, c: {t: 'time'}}));
    expect(items.length, 1);
    expectObject(items[0], {name: 'Biomass', units: 'mg/L', runCount: 1});
  });

  test('dropdown candidates are not unmatched', async () => {
    const full = {path: 't', columns: [{name: 'time'}, {name: 'height'}]};
    const extra = {path: 'u', columns: [{name: 'time'}, {name: 'height'}]};
    const items = unmatched([makeEntry('a', [], [full]), makeEntry('b', [], [full, extra])],
      indexMap({a: {t: 'time'}, b: {t: 'time', u: 'time'}}));
    expect(items.length, 0);
  });

  test('tables without an index are left out', async () => {
    const table = {path: 't', columns: [{name: 'time'}, {name: 'height'}]};
    const items = unmatched([makeEntry('a', [], [table]), makeEntry('b', [], [table])], indexMap({a: {t: 'time'}}));
    expect(items.length, 1);
    expectObject(items[0], {name: 'height', runCount: 1});
  });

  test('values named like a compare row in the same units are left out', async () => {
    const day = (n: number, units = 'g/L') =>
      ({name: 'Final biomass', units, friendlyPath: `Daily cultivation · Day ${n} · Final biomass`});
    const items = unmatched([makeEntry('a', [day(1), day(2), day(3, 'mg/L')]), makeEntry('b', [day(1)])]);
    expect(items.length, 1);
    expectObject(items[0], {name: 'Final biomass', units: 'mg/L', runCount: 1});
  });

  test('one item per name and units, counting its runs', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time'}, {name: 'v', units: 'mL'}]}]),
      makeEntry('b', [], [{path: 't', columns: [{name: 'clock'}, {name: 'V', units: 'ml'}]}]),
    ];
    const items = unmatched(entries, indexMap({a: {t: 'time'}, b: {t: 'clock'}}));
    expect(items.length, 1);
    expectObject(items[0], {name: 'v', units: 'mL', runCount: 2});
  });
});

category('RunComparison: selection gaps', () => {
  test('a scalar missing a charted run is a gap in either selection order', async () => {
    const scalars = matchScalarTargets([
      makeEntry('a', [{name: 'p'}, {name: 'q'}]), makeEntry('b', [{name: 'p'}, {name: 'q'}]),
      makeEntry('c', [{name: 'q'}]),
    ]);
    const p = scalars.find((t) => t.displayName === 'p')!;
    const q = scalars.find((t) => t.displayName === 'q')!;
    expectArray(selectionGaps([p, q]).get(p.key)!, ['c']);
    expectArray(selectionGaps([q, p]).get(p.key)!, ['c']);
    expect(selectionGaps([p, q]).get(q.key)!.length, 0);
  });

  test('a column missing a charted run is a gap in either selection order', async () => {
    const table = (columns: string[]) => ({path: 't', columns: [{name: 'time'}, ...columns.map((name) => ({name}))]});
    const entries = [makeEntry('a', [], [table(['x', 'y'])]), makeEntry('b', [], [table(['x', 'y'])]),
      makeEntry('c', [], [table(['y'])])];
    const targets = matchColumnTargets(entries, indexMap({a: {t: 'time'}, b: {t: 'time'}, c: {t: 'time'}}));
    const x = targets.find((t) => t.displayName === 'x')!;
    const y = targets.find((t) => t.displayName === 'y')!;
    expectArray(selectionGaps([x, y]).get(x.key)!, ['c']);
    expectArray(selectionGaps([y, x]).get(x.key)!, ['c']);
  });

  test('a column picked from another table than the co-target is not a gap', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time'}, {name: 'x'}, {name: 'y'}]}]),
      makeEntry('b', [], [{path: 't', columns: [{name: 'time'}, {name: 'x'}]},
        {path: 'u', columns: [{name: 'time'}, {name: 'y'}]}]),
    ];
    const targets = matchColumnTargets(entries, indexMap({a: {t: 'time'}, b: {t: 'time', u: 'time'}}));
    const x = targets.find((t) => t.displayName === 'x')!;
    const y = targets.find((t) => t.displayName === 'y')!;
    expect(y.bindings.find((b) => b.entryId === 'b')!.tablePath, 'u');
    expect(selectionGaps([x, y]).get(y.key)!.length, 0);
    expect(selectionGaps([y, x]).get(x.key)!.length, 0);
  });
});

category('RunComparison: equal values in mapped units', () => {
  test('converted scalars compare equal in the target units', async () => {
    const convert = buildUnitConverter([{from: 'bar', factor: 100, to: 'kPa'}]);
    const [target] = matchScalarTargets([
      makeEntry('a', [{name: 'P', units: 'bar', value: 1}]), makeEntry('b', [{name: 'P', units: 'kPa', value: 100}]),
    ], undefined, convert);
    expect(isTargetEqualAcrossRuns(target, () => undefined), true);
  });

  test('converted columns compare equal in the target units, indexes compare raw', async () => {
    const convert = buildUnitConverter([{from: 'g/L', factor: 1000, to: 'mg/L'}, {from: 'h', factor: 60, to: 'min'}]);
    const profile = (units: string, timeUnits: string) =>
      ({path: 'profile', columns: [{name: 'Time', units: timeUnits}, {name: 'Biomass', units}]});
    const entries = [makeEntry('a', [], [profile('g/L', 'h')]), makeEntry('b', [], [profile('mg/L', 'min')])];
    const indexes = indexMap({a: {profile: 'Time'}, b: {profile: 'Time'}});
    const [target] = matchColumnTargets(entries, indexes, undefined, undefined, undefined, convert);
    const data: Record<string, Record<string, unknown[]>> = {
      a: {Time: [0, 1], Biomass: [0.1, 2]},
      b: {Time: [0, 1], Biomass: [100, 2000]},
    };
    expect(isTargetEqualAcrossRuns(target, (entryId, _path, column) => data[entryId][column]), true);
    data.b.Time = [0, 60];
    expect(isTargetEqualAcrossRuns(target, (entryId, _path, column) => data[entryId][column]), false);
  });
});

category('RunComparison: units differ status', () => {
  const profile = (column: string, units: string) =>
    ({path: 'profile', columns: [{name: 'Time'}, {name: column, units}]});
  const status = (statuses: ReturnType<typeof getEntryStatuses>, entryId: string) =>
    statuses.find((s) => s.entryId === entryId)!;

  test('a scalar in other units reports units differ with its units', async () => {
    const entries = [
      makeEntry('a', [{name: 'P', units: 'bar'}]), makeEntry('b', [{name: 'P', units: 'bar'}]),
      makeEntry('c', [{name: 'P', units: 'kPa'}]), makeEntry('d', [{name: 'Q', units: 'kPa'}]),
    ];
    const [target] = matchScalarTargets(entries);
    const statuses = getEntryStatuses(entries, target, indexMap({}));
    expect(status(statuses, 'c').reason, 'units differ');
    expect(status(statuses, 'c').units, 'kPa');
    expect(status(statuses, 'd').reason, 'no similar data');
  });

  test('a column in other units reports units differ, unless its index is not set', async () => {
    const entries = [
      makeEntry('a', [], [profile('Biomass', 'g/L')]), makeEntry('b', [], [profile('Biomass', 'g/L')]),
      makeEntry('c', [], [profile('Biomass', 'mg/L')]), makeEntry('d', [], [profile('Biomass', 'mg/L')]),
    ];
    const indexes = indexMap({a: {profile: 'Time'}, b: {profile: 'Time'}, c: {profile: 'Time'}});
    const [target] = matchColumnTargets(entries, indexes);
    const statuses = getEntryStatuses(entries, target, indexes);
    expect(status(statuses, 'c').reason, 'units differ');
    expect(status(statuses, 'c').units, 'mg/L');
    expect(status(statuses, 'd').reason, 'index not set');
  });

  test('units differ honors name mappings and disappears once the units are mapped', async () => {
    const entries = [
      makeEntry('a', [{name: 'P', units: 'bar'}]), makeEntry('b', [{name: 'P', units: 'bar'}]),
      makeEntry('c', [{name: 'pressure', units: 'kPa'}]),
    ];
    const aliases = buildAliasGroups([{from: 'P', to: 'pressure'}]);
    const [target] = matchScalarTargets(entries, aliases);
    expect(status(getEntryStatuses(entries, target, indexMap({})), 'c').reason, 'no similar data');
    expect(status(getEntryStatuses(entries, target, indexMap({}), {aliases}), 'c').reason, 'units differ');
    const convert = buildUnitConverter([{from: 'bar', factor: 100, to: 'kPa'}]);
    const [joined] = matchScalarTargets(entries, aliases, convert);
    expect(status(getEntryStatuses(entries, joined, indexMap({}), {aliases, convert}), 'c').matched, true);
  });
});

category('RunComparison: annotated defaults', () => {
  const annotated = (defaults: Partial<TableNodeInfo>) => {
    const entry = makeEntry('a', [], [{path: 't', columns: [
      {name: 'time', type: 'double'}, {name: 'species', type: 'string'}, {name: 'value'},
    ]}]);
    Object.assign(entry.tables[0], defaults);
    return entry;
  };
  const empty = (): PickerSelections => ({index: {}, split: {}, axis: {}});

  test('empty pickers take the annotated index, split, mode and units', async () => {
    const next = withAnnotatedDefaults(annotated({defaultIndexColumn: 'time', defaultSplitColumn: 'species',
      defaultAxisMode: 'timeseries', defaultTimeUnits: 'h'}), empty());
    expect(next.index.a.t, 'time');
    expect(next.split.a.t, 'species');
    expectObject(next.axis.a.t, {mode: 'timeseries', units: 'h'});
  });

  test('the user\'s picks win', async () => {
    const selections: PickerSelections = {
      index: {a: {t: 'value'}}, split: {a: {t: 'value'}}, axis: {a: {t: {mode: 'points'}}},
    };
    const next = withAnnotatedDefaults(annotated({defaultIndexColumn: 'time', defaultSplitColumn: 'species',
      defaultAxisMode: 'timeseries'}), selections);
    expect(next.index === selections.index && next.split === selections.split && next.axis === selections.axis, true);
  });

  test('defaults that are not valid picks are ignored', async () => {
    const missing = withAnnotatedDefaults(
      annotated({defaultIndexColumn: 'nope', defaultSplitColumn: 'nope'}), empty());
    expect(missing.index.a === undefined && missing.split.a === undefined, true);
    const splitIsIndex = withAnnotatedDefaults(
      annotated({defaultIndexColumn: 'species', defaultSplitColumn: 'species'}), empty());
    expect(splitIsIndex.index.a.t, 'species');
    expect(splitIsIndex.split.a === undefined, true);
    const numericSplit = withAnnotatedDefaults(annotated({defaultSplitColumn: 'value'}), empty());
    expect(numericSplit.split.a === undefined, true);
  });

  test('units are stored for the relative timeseries mode only', async () => {
    const next = withAnnotatedDefaults(annotated({defaultAxisMode: 'points', defaultTimeUnits: 'h'}), empty());
    expect(next.axis.a.t.mode, 'points');
    expect(next.axis.a.t.units === undefined, true);
  });
});
