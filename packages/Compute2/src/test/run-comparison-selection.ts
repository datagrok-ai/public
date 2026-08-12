import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {matchScalarTargets, matchColumnTargets} from '../components/RunComparison/matching';
import {
  getEntryStatuses, matchesFilter, multiValueOverlap, compatibleTargetsFor,
  isSplitCandidate, selectionToMap, computeIndexRows, isTargetEqualAcrossRuns,
  parseTimeUnit, resolveAxisModes, targetAxisMode, resolveDisplayUnit, pickAutoUnit,
} from '../components/RunComparison/selection';
import {ColumnBinding, ColumnTarget, AxisConfigMap, TimeUnit} from '../components/RunComparison/types';
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

  test('multiValueOverlap counts aligned, missing, and conflicting runs', async () => {
    const anchor = columnTarget('anchor', [bind('a', 't'), bind('b', 't'), bind('c', 't2')]);
    const other = columnTarget('other', [bind('a', 't'), bind('c', 't3')]);
    const overlap = multiValueOverlap(anchor, other);
    expect(overlap.aligned, 1);
    expectArray(overlap.missing, ['b']);
    expectArray(overlap.conflicting, ['c']);
  });

  test('zero aligned runs are not suggested', async () => {
    const anchor = columnTarget('anchor', [bind('a', 't')]);
    const other = columnTarget('other', [bind('b', 't')]);
    expect(compatibleTargetsFor(anchor, [anchor, other], () => 'int').length, 1);
  });

  test('unchecking a run keeps the co-target suggested', async () => {
    const table = {path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'height'}, {name: 'velocity'},
    ]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [height] = matchColumnTargets(entries, indexes);
    const edited = matchColumnTargets(entries, indexes, undefined,
      {[height.key]: {'b|t|height': false}});
    const anchor = edited.find((t) => t.displayName === 'height')!;
    expect(anchor.bindings.length, 1);
    expect(compatibleTargetsFor(anchor, edited, () => 'int').length, 2);
  });

  test('cross-table pick becomes conflicting and leaves suggestions', async () => {
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
    expectArray(multiValueOverlap(anchor, velocity).conflicting, ['b']);
    expect(compatibleTargetsFor(anchor, edited, () => 'int').length, 1);
  });

  test('compatibleTargetsFor groups aligned same-table targets', async () => {
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

  test('enabling a sibling swaps the pick and keeps multi-value available', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]}]),
      makeEntry('b', [], [{path: 't', columns: [
        {name: 'time', type: 'int'}, {name: 'height'}, {name: 'heights'},
      ]}]),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [base] = matchColumnTargets(entries, indexes);
    const [anchor] = matchColumnTargets(entries, indexes, undefined,
      {[base.key]: {'b|t|heights': true}});
    expect(anchor.bindings.length, 2);
    expect(anchor.bindings.find((b) => b.entryId === 'b')!.columnName, 'heights');
    expect(compatibleTargetsFor(anchor, [anchor], () => 'int').length, 1);
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

  test('multiValueOverlap on scalars is missing-only, never conflicting', async () => {
    const entries = [
      makeEntry('a', [{name: 'x'}, {name: 'y'}]),
      makeEntry('b', [{name: 'x'}]),
      makeEntry('c', [{name: 'x'}, {name: 'y'}]),
    ];
    const targets = matchScalarTargets(entries);
    const anchor = targets.find((t) => t.displayName === 'x')!;
    const other = targets.find((t) => t.displayName === 'y')!;
    const overlap = multiValueOverlap(anchor, other);
    expect(overlap.aligned, 2);
    expectArray(overlap.missing, ['b']);
    expectArray(overlap.conflicting, []);
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
    expect(parseTimeUnit('m') === undefined, true);
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
    const partial = getEntryStatuses(entries, tgt, indexes, config({a: {t: {units: 's'}}}));
    expect(partial.find((s) => s.entryId === 'a')!.warning === undefined, true);
    expect(partial.find((s) => s.entryId === 'b')!.warning, 'relative timeseries not set');
    const partialScatter = getEntryStatuses(entries, tgt, indexes, config({a: {t: {mode: 'points'}}}));
    expect(partialScatter.find((s) => s.entryId === 'b')!.warning, 'independent points not set');
    const full = getEntryStatuses(entries, tgt, indexes,
      config({a: {t: {units: 's'}}, b: {t: {units: 's'}}}));
    expect(full.every((s) => s.warning === undefined), true);
    const none = getEntryStatuses(entries, tgt, indexes, config({}));
    expect(none.every((s) => s.warning === undefined), true);
  });
});
