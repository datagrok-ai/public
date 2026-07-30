import {category, test, expect, expectFloat, expectArray} from '@datagrok-libraries/test/src/test';
import {
  ComparisonEntryNodes,
  normalizeName, nameSimilarity, nameMatchConfidence, unitsCompatibility, isNumericType,
  matchScalarTargets, matchColumnTargets, getEntryStatuses,
  FUZZY_NAME_THRESHOLD,
} from '../components/RunComparison/comparison-core';

function makeEntry(
  entryId: string,
  scalars: {name: string, valueType?: string, units?: string, value?: number | null}[] = [],
  tables: {path: string, name?: string, columns: {name: string, type?: string, units?: string}[]}[] = [],
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
