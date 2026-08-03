import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {matchScalarTargets, matchColumnTargets} from '../components/RunComparison/matching';
import {
  getEntryStatuses, matchesFilter, multiValueOverlap, compatibleTargetsFor,
  isSplitCandidate, selectionToMap, computeIndexRows,
} from '../components/RunComparison/selection';
import {ColumnBinding, ColumnTarget} from '../components/RunComparison/types';
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
    coverage: bindings.length, defaultCoverage: bindings.length, total: 3,
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
  const anyType = () => true;
  const columns = [{name: 'time', type: 'int'}, {name: 'height'}, {name: 'species', type: 'string'}];

  test('merging groups same-function outputs and intersects columns', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim',
        columns: [...columns, {name: 'extra', type: 'int'}]}]),
    ];
    const rows = computeIndexRows(entries, {}, {}, true, anyType);
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
    const rows = computeIndexRows(entries, {}, {}, true, anyType);
    expect(rows.length, 2);
    expect(rows[0].coverage === undefined, true);
  });

  test('merge off yields one row per table', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
    ];
    expect(computeIndexRows(entries, {}, {}, false, anyType).length, 2);
  });

  test('stale and mixed selections fall back to unset', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
    ];
    const stale = computeIndexRows(entries, {a: {'s/df': 'removed'}, b: {'s/df': 'removed'}}, {}, true, anyType);
    expect(stale[0].current, '');
    const mixed = computeIndexRows(entries, {a: {'s/df': 'time'}, b: {'s/df': 'height'}}, {}, true, anyType);
    expect(mixed[0].current, '');
    const agreed = computeIndexRows(entries, {a: {'s/df': 'time'}, b: {'s/df': 'time'}}, {}, true, anyType);
    expect(agreed[0].current, 'time');
  });

  test('split candidates exclude the chosen index and non-strings', async () => {
    const entries = [
      makeEntry('a', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
      makeEntry('b', [], [{path: 's/df', nqName: 'Pkg:Sim', columns}]),
    ];
    const rows = computeIndexRows(entries, {a: {'s/df': 'time'}, b: {'s/df': 'time'}},
      {a: {'s/df': 'species'}, b: {'s/df': 'species'}}, true, anyType);
    expectArray(rows[0].splitCandidates.map((c) => c.name), ['species']);
    expect(rows[0].currentSplit, 'species');
  });

  test('annotated default index is offered even when its type is toggled off', async () => {
    const gated = (type?: string) => type === 'int';
    const table = {path: 's/df', nqName: 'Pkg:Sim', defaultIndexColumn: 'time_f',
      columns: [{name: 'time_f'}, {name: 'step', type: 'int'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const rows = computeIndexRows(entries, {}, {}, true, gated);
    expectArray(rows[0].candidates.map((c) => c.name), ['time_f', 'step']);
  });

  test('selectionToMap keeps only valid selections', async () => {
    const map = selectionToMap({a: {t: 'time', u: 'bad'}, b: {t: ''}}, (_e, _t, col) => col !== 'bad');
    expect(map.get('a')!.get('t'), 'time');
    expect(map.get('a')!.has('u'), false);
    expect(map.get('b')!.size, 0);
  });
});
