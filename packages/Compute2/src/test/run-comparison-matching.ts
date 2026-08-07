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

category('RunComparison: table compatibility', () => {
  test('index name mismatch blocks clustering', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]}]),
      makeEntry('b', [], [{path: 't', columns: [{name: 'step', type: 'int'}, {name: 'height'}]}]),
    ];
    expect(matchColumnTargets(entries, indexMap({a: {t: 'time'}, b: {t: 'step'}})).length, 0);
  });

  test('split and unsplit tables cluster together', async () => {
    const table = {path: 't', columns: [
      {name: 'time', type: 'int'}, {name: 'height'}, {name: 'species', type: 'string'},
    ]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const targets = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}}),
      indexMap({a: {t: 'species'}}));
    expect(targets.length, 1);
    expect(targets[0].coverage, 2);
    const byEntry = new Map(targets[0].bindings.map((b) => [b.entryId, b.splitColumnName]));
    expect(byEntry.get('a'), 'species');
    expect(byEntry.get('b') == null, true);
  });

  test('differently named splits cluster together', async () => {
    const columns = [
      {name: 'time', type: 'int'}, {name: 'height'},
      {name: 'species', type: 'string'}, {name: 'region', type: 'string'},
    ];
    const entries = [
      makeEntry('a', [], [{path: 't', columns}]),
      makeEntry('b', [], [{path: 't', columns}]),
    ];
    const targets = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}}),
      indexMap({a: {t: 'species'}, b: {t: 'region'}}));
    expect(targets.length, 1);
    expect(targets[0].coverage, 2);
    const byEntry = new Map(targets[0].bindings.map((b) => [b.entryId, b.splitColumnName]));
    expect(byEntry.get('a'), 'species');
    expect(byEntry.get('b'), 'region');
  });

  test('split raw table auto-joins unsplit runs', async () => {
    const runTable = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]};
    const rawTable = {path: 'raw', columns: [
      {name: 'time', type: 'int'}, {name: 'height'}, {name: 'batch', type: 'string'},
    ]};
    const entries = [
      makeEntry('a', [], [runTable]),
      makeEntry('b', [], [runTable]),
      makeEntry('r', [], [rawTable], 'raw'),
    ];
    const targets = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}, r: {raw: 'time'}}),
      indexMap({r: {raw: 'batch'}}));
    expect(targets.length, 1);
    expect(targets[0].coverage, 3);
    const rawCandidate = targets[0].candidates.find((c) => c.binding.entryId === 'r')!;
    expect(rawCandidate.enabled, true);
    expect(rawCandidate.binding.splitColumnName, 'batch');
  });
});

category('RunComparison: candidates', () => {
  test('other compatible columns attach as disabled candidates', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'temperature'}]}]),
      makeEntry('b', [], [{path: 't', columns: [
        {name: 'time', type: 'int'}, {name: 'temperature'}, {name: 'temperatures'},
      ]}]),
    ];
    const targets = matchColumnTargets(entries, indexMap({a: {t: 'time'}, b: {t: 'time'}}));
    expect(targets.length, 1);
    const [target] = targets;
    expect(target.candidates.length, 3);
    const extra = target.candidates.find((c) => c.binding.columnName === 'temperatures')!;
    expect(extra.auto, false);
    expect(extra.enabled, false);
    expect(extra.confidence, 'fuzzy');
    expect(target.candidates.filter((c) => c.auto).every((c) => c.enabled), true);
    expect(target.bindings.length, 2);
    expect(target.coverage, 2);
    expect(target.confidence, 'exact');
  });

  test('raw items are enabled in every compatible cluster', async () => {
    const tables = [
      {path: 's1/df', name: 'Heating', columns: [{name: 'time', type: 'int'}, {name: 'height'}]},
      {path: 's2/df', name: 'Cooling', columns: [{name: 'time', type: 'int'}, {name: 'height'}]},
    ];
    const entries = [
      makeEntry('a', [], tables),
      makeEntry('b', [], tables),
      makeEntry('raw', [], [{path: 'exp', columns: [{name: 'time', type: 'int'}, {name: 'height'}]}], 'raw'),
    ];
    const targets = matchColumnTargets(entries, indexMap({
      a: {'s1/df': 'time', 's2/df': 'time'},
      b: {'s1/df': 'time', 's2/df': 'time'},
      raw: {exp: 'time'},
    }));
    expect(targets.length, 2);
    for (const target of targets) {
      const rawCandidate = target.candidates.find((c) => c.binding.entryId === 'raw')!;
      expect(rawCandidate.enabled, true);
      expect(target.coverage, 3);
    }
    expect(targets.filter((t) =>
      t.candidates.find((c) => c.binding.entryId === 'raw')!.auto).length, 1);
  });

  test('raw attachment resurrects a single-entry cluster', async () => {
    const entries = [
      makeEntry('a', [], [
        {path: 't1', name: 'T1', columns: [{name: 'time', type: 'int'}, {name: 'height'}]},
        {path: 't2', name: 'T2', columns: [{name: 'time', type: 'int'}, {name: 'height'}]},
      ]),
      makeEntry('raw', [], [{path: 'exp', columns: [{name: 'time', type: 'int'}, {name: 'height'}]}], 'raw'),
    ];
    const targets = matchColumnTargets(entries,
      indexMap({a: {t1: 'time', t2: 'time'}, raw: {exp: 'time'}}));
    expect(targets.length, 2);
    for (const target of targets)
      expect(target.coverage, 2);
  });

  test('raw items enable only on non-fuzzy matches', async () => {
    const wfTable = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'temperature'}]};
    const entries = [
      makeEntry('a', [], [wfTable]),
      makeEntry('b', [], [wfTable]),
      makeEntry('raw', [], [{path: 'exp', columns: [
        {name: 'time', type: 'int'}, {name: 'Temperature'}, {name: 'temperatures'},
      ]}], 'raw'),
    ];
    const targets = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}, raw: {exp: 'time'}}));
    const target = targets.find((t) => t.displayName === 'temperature')!;
    const normalized = target.candidates.find((c) => c.binding.columnName === 'Temperature')!;
    expect(normalized.enabled, true);
    const fuzzy = target.candidates.find((c) => c.binding.columnName === 'temperatures')!;
    expect(fuzzy.enabled, false);
    expect(target.bindings.filter((b) => b.entryId === 'raw').length, 1);
  });

  test('fuzzy greedy raw pick stays unchecked', async () => {
    const wfTable = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'temperature'}]};
    const entries = [
      makeEntry('a', [], [wfTable]),
      makeEntry('b', [], [wfTable]),
      makeEntry('raw', [], [{path: 'exp', columns: [
        {name: 'time', type: 'int'}, {name: 'temperatures'},
      ]}], 'raw'),
    ];
    const [target] = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}, raw: {exp: 'time'}}));
    const rawCandidate = target.candidates.find((c) => c.binding.entryId === 'raw')!;
    expect(rawCandidate.auto, true);
    expect(rawCandidate.enabled, false);
    expect(target.coverage, 2);
  });

  test('overrides flip candidates and update derived fields', async () => {
    const entries = [
      makeEntry('a', [], [{path: 't', columns: [{name: 'time', type: 'int'}, {name: 'temperature'}]}]),
      makeEntry('b', [], [{path: 't', columns: [
        {name: 'time', type: 'int'}, {name: 'temperature'}, {name: 'temperatures', units: 'K'},
      ]}]),
    ];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [base] = matchColumnTargets(entries, indexes);
    const key = base.key;

    const [disabled] = matchColumnTargets(entries, indexes, undefined, {[key]: {'b|t|temperature': false}});
    expect(disabled.bindings.length, 1);
    expect(disabled.coverage, 1);
    expect(disabled.confidence, 'exact');

    // radio semantics: the explicit pick replaces the run's auto one
    const [swapped] = matchColumnTargets(entries, indexes, undefined, {[key]: {'b|t|temperatures': true}});
    expect(swapped.bindings.length, 2);
    expect(swapped.coverage, 2);
    expect(swapped.confidence, 'fuzzy');
    expect(swapped.unitsWarning, true);
    expect(swapped.candidates.find((c) => c.binding.columnName === 'temperature' &&
      c.binding.entryId === 'b')!.enabled, false);

    const [stale] = matchColumnTargets(entries, indexes, undefined, {[key]: {'zzz|t|nope': false}});
    expect(stale.bindings.length, 2);
  });

  test('at most one candidate per run is enabled by default', async () => {
    const wfTable = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'temperature'}]};
    const entries = [
      makeEntry('a', [], [wfTable]),
      makeEntry('b', [], [wfTable]),
      makeEntry('raw', [], [{path: 'exp', columns: [
        {name: 'time', type: 'int'}, {name: 'Temperature'}, {name: 'TEMPERATURE'},
      ]}], 'raw'),
    ];
    const [target] = matchColumnTargets(entries,
      indexMap({a: {t: 'time'}, b: {t: 'time'}, raw: {exp: 'time'}}));
    const rawCandidates = target.candidates.filter((c) => c.binding.entryId === 'raw');
    expect(rawCandidates.length, 2);
    expect(rawCandidates.filter((c) => c.enabled).length, 1);
  });

  test('user toggles never remove the target', async () => {
    const table = {path: 't', columns: [{name: 'time', type: 'int'}, {name: 'height'}]};
    const entries = [makeEntry('a', [], [table]), makeEntry('b', [], [table])];
    const indexes = indexMap({a: {t: 'time'}, b: {t: 'time'}});
    const [base] = matchColumnTargets(entries, indexes);
    const [target] = matchColumnTargets(entries, indexes, undefined,
      {[base.key]: {'a|t|height': false, 'b|t|height': false}});
    expect(target.key, base.key);
    expect(target.bindings.length, 0);
    expect(target.coverage, 0);
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
