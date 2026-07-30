import {category, test, expect, expectFloat, expectArray} from '@datagrok-libraries/test/src/test';
import {
  ComparisonEntryNodes,
  normalizeName, nameSimilarity, nameMatchConfidence, unitsCompatibility, isNumericType,
  matchScalarTargets, matchColumnTargets, getEntryStatuses,
  matchesFilter, bindingSignature, compatibleTargetsFor, isSplitCandidate, selectionToMap, computeIndexRows,
  FUZZY_NAME_THRESHOLD,
} from '../components/RunComparison/comparison-core';

function makeEntry(
  entryId: string,
  scalars: {name: string, valueType?: string, units?: string, value?: number | null}[] = [],
  tables: {
    path: string, name?: string, nqName?: string,
    columns: {name: string, type?: string, units?: string}[],
  }[] = [],
): ComparisonEntryNodes {
  return {
    entryId,
    entryName: entryId,
    scalars: scalars.map((s, i) => ({
      path: `io${i}`,
      name: s.name,
      valueType: s.valueType ?? 'double',
      units: s.units,
      value: s.value ?? 0,
    })),
    tables: tables.map((t) => ({
      path: t.path,
      name: t.name ?? t.path,
      nqName: t.nqName,
      columns: t.columns.map((c) => ({name: c.name, type: c.type ?? 'double', units: c.units})),
      rowCount: 10,
    })),
  };
}

const indexMap = (spec: Record<string, Record<string, string>>) =>
  new Map(Object.entries(spec).map(([entryId, tables]) => [entryId, new Map(Object.entries(tables))]));

category('RunComparison: names and units', () => {
  test('normalizeName collapses case and separators', async () => {
    expect(normalizeName('  Init_Temp '), 'init temp');
    expect(normalizeName('init-temp'), 'init temp');
    expect(normalizeName('Init  Temp'), 'init temp');
  });

  test('nameSimilarity is 1 for equal normalized names', async () => {
    expectFloat(nameSimilarity('Init_Temp', 'init temp'), 1);
  });

  test('nameSimilarity is high for near names and low for unrelated', async () => {
    expect(nameSimilarity('concentration', 'concentrations') >= FUZZY_NAME_THRESHOLD, true);
    expect(nameSimilarity('temperature', 'pressure') < FUZZY_NAME_THRESHOLD, true);
  });

  test('nameMatchConfidence tiers', async () => {
    expect(nameMatchConfidence('temp', 'temp'), 'exact');
    expect(nameMatchConfidence('Init_Temp', 'init temp'), 'normalized');
    expect(nameMatchConfidence('concentration', 'concentrations'), 'fuzzy');
    expect(nameMatchConfidence('temperature', 'pressure'), null);
  });

  test('unitsCompatibility match, warn, mismatch', async () => {
    expect(unitsCompatibility('mg', 'MG'), 'match');
    expect(unitsCompatibility(undefined, undefined), 'match');
    expect(unitsCompatibility('mg', undefined), 'warn');
    expect(unitsCompatibility('mg', 'g'), 'mismatch');
  });

  test('isNumericType accepts numbers only', async () => {
    expect(isNumericType('int'), true);
    expect(isNumericType('double'), true);
    expect(isNumericType('string'), false);
    expect(isNumericType('datetime'), false);
  });
});

category('RunComparison: scalar matching', () => {
  test('groups scalars with normalized names across entries', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'Init_Temp', value: 1}]),
      makeEntry('b', [{name: 'init temp', value: 2}]),
    ]);
    expect(targets.length, 1);
    expect(targets[0].confidence, 'normalized');
    expect(targets[0].coverage, 2);
    expectArray(targets[0].bindings.map((b) => b.value), [1, 2]);
  });

  test('ignores non-numeric scalars', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'label', valueType: 'string'}]),
      makeEntry('b', [{name: 'label', valueType: 'string'}]),
    ]);
    expect(targets.length, 0);
  });

  test('does not group scalars with mismatching units', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'mass', units: 'mg'}]),
      makeEntry('b', [{name: 'mass', units: 'g'}]),
    ]);
    expect(targets.length, 0);
  });

  test('warns when units are partially missing', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'mass', units: 'mg'}]),
      makeEntry('b', [{name: 'mass'}]),
    ]);
    expect(targets.length, 1);
    expect(targets[0].unitsWarning, true);
  });

  test('drops groups present in a single entry only', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'unique'}, {name: 'shared'}]),
      makeEntry('b', [{name: 'shared'}]),
    ]);
    expect(targets.length, 1);
    expect(targets[0].displayName, 'shared');
    expect(targets[0].confidence, 'exact');
  });

  test('takes at most one scalar per entry into a target', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'temp', value: 1}, {name: 'temp2', value: 10}]),
      makeEntry('b', [{name: 'temp', value: 2}]),
    ]);
    const shared = targets.find((t) => t.displayName === 'temp')!;
    expect(shared.bindings.length, 2);
    expectArray(shared.bindings.map((b) => b.value), [1, 2]);
  });
});

category('RunComparison: column matching', () => {
  const tableA = {path: 'sim/out', name: 'Cooling', columns: [
    {name: 'time'}, {name: 'temperature'}, {name: 'label', type: 'string'},
  ]};
  const tableB = {path: 'result', name: 'Cooling result', columns: [
    {name: 'time'}, {name: 'temperatures'},
  ]};

  test('requires a user-defined index on both tables', async () => {
    const entries = [makeEntry('a', [], [tableA]), makeEntry('b', [], [tableB])];
    expect(matchColumnTargets(entries, indexMap({})).length, 0);
    expect(matchColumnTargets(entries, indexMap({a: {'sim/out': 'time'}})).length, 0);
    const targets = matchColumnTargets(entries, indexMap({a: {'sim/out': 'time'}, b: {result: 'time'}}));
    expect(targets.length, 1);
    expect(targets[0].displayName, 'temperature');
    expect(targets[0].confidence, 'fuzzy');
    expect(targets[0].coverage, 2);
  });

  test('excludes the index column and non-numeric columns from candidates', async () => {
    const entries = [makeEntry('a', [], [tableA]), makeEntry('b', [], [tableB])];
    const targets = matchColumnTargets(entries, indexMap({a: {'sim/out': 'time'}, b: {result: 'time'}}));
    const names = targets.map((t) => t.displayName);
    expect(names.includes('time'), false);
    expect(names.includes('label'), false);
  });

  test('bindings carry table path and index column', async () => {
    const entries = [makeEntry('a', [], [tableA]), makeEntry('b', [], [tableB])];
    const [target] = matchColumnTargets(entries, indexMap({a: {'sim/out': 'time'}, b: {result: 'time'}}));
    const bindingA = target.bindings.find((b) => b.entryId === 'a')!;
    expect(bindingA.tablePath, 'sim/out');
    expect(bindingA.indexColumnName, 'time');
    expect(bindingA.columnName, 'temperature');
  });

  test('same-named clusters get unique keys', async () => {
    const columns = [{name: 'time', type: 'int'}, {name: 'height'}];
    const tables = [
      {path: 'step1/df', name: 'Results', columns},
      {path: 'step2/df', name: 'Results', columns},
    ];
    const entries = [makeEntry('a', [], tables), makeEntry('b', [], tables)];
    const index = indexMap({
      a: {'step1/df': 'time', 'step2/df': 'time'},
      b: {'step1/df': 'time', 'step2/df': 'time'},
    });
    const targets = matchColumnTargets(entries, index);
    expect(targets.length, 2);
    expect(new Set(targets.map((t) => t.key)).size, 2);
  });
});

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
});

category('RunComparison: split matching', () => {
  const table = {path: 't', columns: [
    {name: 'time', type: 'int'}, {name: 'height'}, {name: 'species', type: 'string'},
  ]};
  const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
  const splits = indexMap({a: {t: 'species'}, b: {t: 'species'}});

  test('split column is excluded from candidates and carried into bindings', async () => {
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const targets = matchColumnTargets(entries, indexes, splits);
    expect(targets.length, 1);
    expect(targets[0].displayName, 'height');
    expect(targets[0].bindings.every((b) => b.splitColumnName === 'species'), true);
  });

  test('split choice changes the bindings signature', async () => {
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const [withSplit] = matchColumnTargets(entries, indexes, splits);
    const [withoutSplit] = matchColumnTargets(entries, indexes);
    expect(bindingSignature(withSplit) === bindingSignature(withoutSplit), false);
  });

  test('isSplitCandidate requires a string column that is not the index', async () => {
    expect(isSplitCandidate({name: 'species', type: 'string'}, 'time'), true);
    expect(isSplitCandidate({name: 'time', type: 'string'}, 'time'), false);
    expect(isSplitCandidate({name: 'height', type: 'double'}, 'time'), false);
  });
});

category('RunComparison: cluster pairing', () => {
  test('same-named columns pair by table name across entries', async () => {
    const tables = (suffix: string) => [
      {path: `s1/df${suffix}`, name: 'Heating', columns: [{name: 'time', type: 'int'}, {name: 'height'}]},
      {path: `s2/df${suffix}`, name: 'Cooling', columns: [{name: 'time', type: 'int'}, {name: 'height'}]},
    ];
    const entries = [makeEntry('a', [], tables('A')), makeEntry('b', [], tables('B'))];
    const indexes = indexMap({
      a: {'s1/dfA': 'time', 's2/dfA': 'time'},
      b: {'s1/dfB': 'time', 's2/dfB': 'time'},
    });
    const targets = matchColumnTargets(entries, indexes);
    expect(targets.length, 2);
    for (const target of targets) {
      const names = new Set(target.bindings.map((b) => b.tableName));
      expect(names.size, 1);
    }
  });
});

category('RunComparison: entry statuses extra', () => {
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

  test('compatibleTargetsFor groups by full signature', async () => {
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

  test('selectionToMap keeps only valid selections', async () => {
    const map = selectionToMap({a: {t: 'time', u: 'bad'}, b: {t: ''}}, (_e, _t, col) => col !== 'bad');
    expect(map.get('a')!.get('t'), 'time');
    expect(map.get('a')!.has('u'), false);
    expect(map.get('b')!.size, 0);
  });
});
