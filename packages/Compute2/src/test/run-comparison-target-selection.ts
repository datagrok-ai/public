import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {
  TargetSelection, EMPTY_TARGET_SELECTION, isCharted, clickTarget, setMultiMode, reconcileTargets, pruneOverrides,
  chartedTargets, withCandidateEnabled, withoutTargetOverrides,
} from '../components/RunComparison/target-selection';
import {ColumnTarget, ScalarTarget, ComparisonTarget, ColumnBinding} from '../components/RunComparison/types';

// a and b chart together, c has no partners
const partners: Record<string, string[]> = {a: ['a', 'b'], b: ['a', 'b'], c: ['c']};
const compatible = (key: string) => partners[key] ?? [];

const single = (selectedKey: string): TargetSelection => ({selectedKey, multiMode: false, multiKeys: []});
const multi = (selectedKey: string, multiKeys: string[]): TargetSelection =>
  ({selectedKey, multiMode: true, multiKeys});

const expectState = (actual: TargetSelection, expected: TargetSelection) => {
  expect(actual.selectedKey, expected.selectedKey);
  expect(actual.multiMode, expected.multiMode);
  expectArray(actual.multiKeys, expected.multiKeys);
};

category('RunComparison: target selection clicks', () => {
  test('click selects, a second click clears', async () => {
    const selected = clickTarget(EMPTY_TARGET_SELECTION, 'a', false, compatible);
    expectState(selected, single('a'));
    expectState(clickTarget(selected, 'a', false, compatible), single(''));
  });

  test('click on another value replaces the selection', async () => {
    expectState(clickTarget(single('a'), 'c', false, compatible), single('c'));
  });

  test('shift click on a compatible value enters multi mode with both', async () => {
    expectState(clickTarget(single('a'), 'b', true, compatible), multi('a', ['a', 'b']));
  });

  test('shift click on an incompatible value with partners re-anchors in multi mode', async () => {
    const state = clickTarget(single('c'), 'a', true, compatible);
    expectState(state, multi('a', ['a']));
  });

  test('shift click on a value without partners only selects it', async () => {
    expectState(clickTarget(single('a'), 'c', true, compatible), single('c'));
  });

  test('shift click without a selection anchors on the clicked value', async () => {
    expectState(clickTarget(EMPTY_TARGET_SELECTION, 'a', true, compatible), multi('a', ['a']));
    expectState(clickTarget(EMPTY_TARGET_SELECTION, 'c', true, compatible), single('c'));
  });

  test('shift click on the selected value enters multi mode with it alone', async () => {
    expectState(clickTarget(single('a'), 'a', true, compatible), multi('a', ['a']));
  });

  test('click in multi mode adds and removes values', async () => {
    const added = clickTarget(multi('a', ['a']), 'b', false, compatible);
    expectState(added, multi('a', ['a', 'b']));
    expectState(clickTarget(added, 'b', true, compatible), multi('a', ['a']));
  });

  test('unchecking the anchor moves it to the first charted value', async () => {
    expectState(clickTarget(multi('a', ['a', 'b', 'c']), 'a', false, compatible), multi('b', ['b', 'c']));
  });

  test('unchecking the last value exits multi mode with it selected', async () => {
    expectState(clickTarget(multi('a', ['a']), 'a', false, compatible), single('a'));
  });

  test('charted values follow the mode', async () => {
    expect(isCharted(single('a'), 'a'), true);
    expect(isCharted(single('a'), 'b'), false);
    expect(isCharted(multi('a', ['a', 'b']), 'b'), true);
    expect(isCharted(multi('a', ['a', 'b']), 'c'), false);
    expect(isCharted({selectedKey: 'a', multiMode: false, multiKeys: ['b']}, 'b'), false);
  });
});

category('RunComparison: target selection mode toggle', () => {
  test('entering multi mode seeds it with the selected value', async () => {
    expectState(setMultiMode(single('a'), true), multi('a', ['a']));
    expectState(setMultiMode(EMPTY_TARGET_SELECTION, true), multi('', []));
  });

  test('leaving multi mode keeps the first charted value selected', async () => {
    const state = setMultiMode(multi('b', ['a', 'b']), false);
    expect(state.selectedKey, 'a');
    expect(state.multiMode, false);
    expect(isCharted(state, 'a'), true);
    expect(isCharted(state, 'b'), false);
  });

  test('leaving an empty multi mode keeps the selection', async () => {
    expectState(setMultiMode(multi('a', []), false), {selectedKey: 'a', multiMode: false, multiKeys: []});
  });

  test('a toggle echo returns the same state', async () => {
    const state = multi('a', ['a', 'b']);
    expect(setMultiMode(state, true) === state, true);
    const plain = single('a');
    expect(setMultiMode(plain, false) === plain, true);
  });
});

category('RunComparison: target selection reconciliation', () => {
  test('a vanished value is dropped and the anchor moves to a charted one', async () => {
    expectState(reconcileTargets(multi('a', ['a', 'b', 'c']), new Set(['b', 'c'])), multi('b', ['b', 'c']));
  });

  test('a vanished non-anchor value keeps the anchor', async () => {
    expectState(reconcileTargets(multi('a', ['a', 'b']), new Set(['a'])), multi('a', ['a']));
  });

  test('losing every value exits multi mode', async () => {
    const state = reconcileTargets(multi('a', ['a', 'b']), new Set(['c']));
    expect(state.multiMode, false);
    expectArray(state.multiKeys, []);
  });

  test('a vanished single selection is kept for when the value returns', async () => {
    const state = single('a');
    expect(reconcileTargets(state, new Set(['b'])) === state, true);
  });

  test('unchanged targets return the same state', async () => {
    const state = multi('a', ['a', 'b']);
    expect(reconcileTargets(state, new Set(['a', 'b', 'c'])) === state, true);
  });
});

const binding = (entryId: string, columnName: string): ColumnBinding => ({
  entryId, tablePath: 't', tableName: 't', columnName, indexColumnName: 'time',
});

const columnTarget = (key: string, bindings: ColumnBinding[]): ColumnTarget => ({
  kind: 'column', key, displayName: key, confidence: 'exact', unitsWarning: false,
  coverage: bindings.length, total: bindings.length, bindings,
  candidates: bindings.map((b) => ({binding: b, confidence: 'exact', unitsWarn: false, auto: true, enabled: true})),
});

const scalarTarget = (key: string): ScalarTarget => ({
  kind: 'scalar', key, displayName: key, confidence: 'exact', unitsWarning: false,
  coverage: 2, total: 2, bindings: [],
});

category('RunComparison: candidate override pruning', () => {
  const targets: ComparisonTarget[] = [
    columnTarget('column:x', [binding('r1', 'x'), binding('r2', 'x'), binding('r2', 'x2')]),
    scalarTarget('scalar:y'),
  ];

  test('stale candidates and targets are dropped, valid ones kept', async () => {
    const pruned = pruneOverrides({
      'column:x': {'r1|t|x': false, 'r3|t|x': true},
      'column:gone': {'r1|t|x': true},
      'scalar:y': {'r1|t|y': true},
    }, targets);
    expectArray(Object.keys(pruned), ['column:x']);
    expectArray(Object.keys(pruned['column:x']), ['r1|t|x']);
    expect(pruned['column:x']['r1|t|x'], false);
  });

  test('a target left without valid overrides is removed', async () => {
    expectArray(Object.keys(pruneOverrides({'column:x': {'r3|t|x': true}}, targets)), []);
  });

  test('valid overrides return the same object', async () => {
    const overrides = {'column:x': {'r1|t|x': false, 'r2|t|x2': true}};
    expect(pruneOverrides(overrides, targets) === overrides, true);
  });
});

category('RunComparison: charted targets', () => {
  const targets: ComparisonTarget[] = [
    scalarTarget('s1'), scalarTarget('s2'), scalarTarget('s3'), columnTarget('c1', [binding('r1', 'x')]),
  ];
  const keys = (list: ComparisonTarget[]) => list.map((target) => target.key);

  test('single mode charts the selection only', async () => {
    expectArray(keys(chartedTargets(single('s2'), targets)), ['s2']);
    expectArray(keys(chartedTargets(single(''), targets)), []);
    expectArray(keys(chartedTargets(single('gone'), targets)), []);
  });

  test('multi mode charts the anchor first, then list order', async () => {
    expectArray(keys(chartedTargets(multi('s3', ['s1', 's3', 's2']), targets)), ['s3', 's1', 's2']);
  });

  test('targets of another kind than the anchor are left out', async () => {
    expectArray(keys(chartedTargets(multi('s1', ['s1', 'c1']), targets)), ['s1']);
  });

  test('a vanished anchor charts nothing', async () => {
    expectArray(keys(chartedTargets(multi('gone', ['s1']), targets)), []);
  });
});

category('RunComparison: candidate toggles', () => {
  const target = columnTarget('column:x', [binding('r1', 'x'), binding('r2', 'x'), binding('r2', 'x2')]);
  target.candidates[2].enabled = false;
  const [r1x, r2x, r2x2] = target.candidates;

  test('enabling a candidate disables the run\'s enabled one only', async () => {
    const next = withCandidateEnabled({}, target, r2x2, true);
    expectArray(Object.keys(next['column:x']).sort(), ['r2|t|x', 'r2|t|x2']);
    expect(next['column:x']['r2|t|x2'], true);
    expect(next['column:x']['r2|t|x'], false);
  });

  test('disabling a candidate leaves the other runs alone', async () => {
    const next = withCandidateEnabled({}, target, r1x, false);
    expectArray(Object.keys(next['column:x']), ['r1|t|x']);
    expect(next['column:x']['r1|t|x'], false);
  });

  test('toggles keep the overrides of other targets and earlier toggles', async () => {
    const next = withCandidateEnabled({'column:y': {'r1|t|y': false}, 'column:x': {'r1|t|x': false}},
      target, r2x, false);
    expect(next['column:y']['r1|t|y'], false);
    expect(next['column:x']['r1|t|x'], false);
    expect(next['column:x']['r2|t|x'], false);
  });

  test('reset drops the target\'s overrides only', async () => {
    const next = withoutTargetOverrides({'column:x': {'r1|t|x': false}, 'column:y': {'r1|t|y': false}}, 'column:x');
    expectArray(Object.keys(next), ['column:y']);
  });
});
