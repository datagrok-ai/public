import {category, test, expect, expectFloat, expectArray} from '@datagrok-libraries/test/src/test';
import {isNumericType} from '../components/RunComparison/types';
import {
  normalizeName, nameSimilarity, nameMatchConfidence, unitsCompatibility,
  matchScalarTargets, matchColumnTargets, FUZZY_NAME_THRESHOLD,
} from '../components/RunComparison/matching';
import {makeEntry, indexMap} from './run-comparison-fixtures';

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

  test('does not group columns with mismatching units', async () => {
    const table = (units: string) => ({path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'mass', units},
    ]});
    const entries = [makeEntry('a', [], [table('mg')]), makeEntry('b', [], [table('g')])];
    const targets = matchColumnTargets(entries, indexMap({a: {t: 'time'}, b: {t: 'time'}}));
    expect(targets.length, 0);
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

  test('split column is excluded from candidates and carried into bindings', async () => {
    const table = {path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'height'}, {name: 'species', type: 'string'},
    ]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const targets = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}}),
      indexMap({a: {t: 'species'}, b: {t: 'species'}}));
    expect(targets.length, 1);
    expect(targets[0].displayName, 'height');
    expect(targets[0].bindings.every((b) => b.splitColumnName === 'species'), true);
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

  test('exact-name cluster wins over an earlier fuzzy cluster', async () => {
    const targets = matchScalarTargets([
      makeEntry('a', [{name: 'temperatures', value: 1}, {name: 'temperature', value: 2}]),
      makeEntry('b', [{name: 'temperature', value: 3}]),
    ]);
    expect(targets.length, 1);
    expect(targets[0].displayName, 'temperature');
    expect(targets[0].confidence, 'exact');
    expectArray(targets[0].bindings.map((b) => b.value), [2, 3]);
  });
});
