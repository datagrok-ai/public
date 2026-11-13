import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {
  ALL_INPUTS_TOGGLES, DEFAULT_INPUT_TOGGLES, getInputToggles, InputToggle,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {RestrictionType} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {getToggledInputs} from '../components/RFV/input-toggles';

const prop = (name: string, category?: string) => DG.Property.js(name, DG.TYPE.FLOAT, {category});
const consistency = (restriction: RestrictionType, inconsistent = false): ConsistencyInfo =>
  ({restriction, inconsistent, assignedValue: 1});

const PROPS = ['disabled', 'restricted', 'info', 'none', 'untracked'].map((name) => prop(name));
const RESTRICTIONS: Record<string, ConsistencyInfo> = {
  disabled: consistency('disabled'),
  restricted: consistency('restricted'),
  info: consistency('info'),
  none: consistency('none'),
};

const toggled = (
  props: DG.Property[], toggles: InputToggle[], shown: string[] = [], hiddenByMeta: string[] = [],
  validations: Record<string, any> = {}, consistencies = RESTRICTIONS,
) => {
  const {hidden, counts} =
    getToggledInputs(props, toggles, new Set(shown), new Set(hiddenByMeta), validations, consistencies);
  return {hidden: [...hidden].sort(), counts: Object.fromEntries(counts)};
};

category('RFV: input toggles', () => {
  test('The default toggle hides disabled and restricted inputs only', async () => {
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES),
      {hidden: ['disabled', 'restricted'], counts: {'linked inputs': 2}});
  });

  test('A shown toggle keeps its count', async () => {
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES, ['linked inputs']),
      {hidden: [], counts: {'linked inputs': 2}});
  });

  test('Inputs hidden by meta are left out', async () => {
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES, [], ['disabled']),
      {hidden: ['restricted'], counts: {'linked inputs': 1}});
  });

  test('The keepVisible annotation keeps an input visible', async () => {
    const script = DG.Script.create(['//name: t', '//language: javascript',
      '//input: double x {keepVisible: true}', '//input: double y'].join('\n'));
    expectDeepEqual(toggled(script.inputs, DEFAULT_INPUT_TOGGLES, [], [], {},
      {x: consistency('disabled'), y: consistency('disabled')}),
    {hidden: ['y'], counts: {'linked inputs': 1}});
  });

  test('Errors, warnings and inconsistent values reveal an input', async () => {
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES, [], [], {disabled: {errors: ['bad']}}),
      {hidden: ['restricted'], counts: {'linked inputs': 1}});
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES, [], [], {restricted: {warnings: ['careful']}}),
      {hidden: ['disabled'], counts: {'linked inputs': 1}});
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES, [], [], {restricted: {notifications: ['note']}}),
      {hidden: ['disabled', 'restricted'], counts: {'linked inputs': 2}});
    expectDeepEqual(toggled(PROPS, DEFAULT_INPUT_TOGGLES, [], [], {},
      {...RESTRICTIONS, restricted: consistency('restricted', true)}),
    {hidden: ['disabled'], counts: {'linked inputs': 1}});
  });

  test('Toggles select by category and input name, the first match wins', async () => {
    const props = [prop('pump', 'Pump'), prop('valve', 'Valves'), prop('atol', 'Solver'), prop('rtol'), prop('other')];
    const toggles: InputToggle[] = [
      {caption: 'device params', categories: ['Pump', 'Valves'], linked: true},
      {caption: 'solver settings', inputs: ['atol', 'rtol']},
    ];
    const consistencies = {pump: consistency('restricted'), valve: consistency('info'), other: consistency('disabled')};
    expectDeepEqual(toggled(props, toggles, ['solver settings'], [], {}, consistencies),
      {hidden: ['pump'], counts: {'device params': 1, 'solver settings': 2}});
  });

  test('An entry without selectors hides every input but those with errors', async () => {
    expectDeepEqual(toggled(PROPS, [{caption: 'inputs'}], [], [], {info: {errors: ['bad']}}),
      {hidden: ['disabled', 'none', 'restricted', 'untracked'], counts: {inputs: 4}});
  });

  test('inputsHidden puts every input behind the linked inputs row, inputToggles wins', async () => {
    const script = (...meta: string[]) =>
      DG.Script.create(['//name: t', '//language: javascript', ...meta, '//input: double x'].join('\n'));
    expectDeepEqual(getInputToggles(script()), DEFAULT_INPUT_TOGGLES);
    expectDeepEqual(getInputToggles(script('//meta.inputsHidden: true')), ALL_INPUTS_TOGGLES);
    expectDeepEqual(getInputToggles(script('//meta.inputsHidden: true', '//meta.inputToggles: []')), []);
    expectDeepEqual(toggled(PROPS, ALL_INPUTS_TOGGLES),
      {hidden: ['disabled', 'info', 'none', 'restricted', 'untracked'], counts: {'linked inputs': 5}});
  });

  test('Malformed or non-array inputToggles falls back to the default', async () => {
    const script = (meta: string) =>
      DG.Script.create(['//name: t', '//language: javascript', meta, '//input: double x'].join('\n'));
    expectDeepEqual(getInputToggles(script('//meta.inputToggles: [{caption: x}]')), DEFAULT_INPUT_TOGGLES);
    expectDeepEqual(getInputToggles(script('//meta.inputToggles: {"caption": "x"}')), DEFAULT_INPUT_TOGGLES);
  });

  test('Entries with the same caption share a row', async () => {
    const toggles = [...DEFAULT_INPUT_TOGGLES, {caption: 'linked inputs', inputs: ['untracked']}];
    expectDeepEqual(toggled(PROPS, toggles),
      {hidden: ['disabled', 'restricted', 'untracked'], counts: {'linked inputs': 3}});
  });

  test('An empty list turns toggles off', async () => {
    expectDeepEqual(toggled(PROPS, []), {hidden: [], counts: {}});
  });
});
