import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BehaviorSubject} from 'rxjs';
import {category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {InputForm} from '@datagrok-libraries/webcomponents-vue';
import {applyInputsVisibility} from '@datagrok-libraries/webcomponents-vue/src/InputForm/utils';
import {
  ALL_INPUTS_TOGGLES, DEFAULT_INPUT_TOGGLES, getInputToggles, InputToggle,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {RestrictionType} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {getToggledInputs} from '../components/RFV/input-toggles';
import {awaitWebComponents} from './utils';

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

class FakeInput {
  writes: boolean[] = [];
  root = undefined as unknown as HTMLElement;
  private shown = true;
  get visible() {
    return this.shown;
  }
  set visible(value: boolean) {
    this.shown = value;
    this.writes.push(value);
  }
}

category('RFV: input visibility', () => {
  test('Visibility is written only when it changes', async () => {
    const inputs: Record<string, FakeInput> = {a: new FakeInput(), b: new FakeInput(), c: new FakeInput()};
    const apply = (hidden: string[], meta: Record<string, any> = {}) =>
      applyInputsVisibility([...Object.keys(inputs), 'missing'], (name) => inputs[name], meta, new Set(hidden));
    const writes = () => Object.fromEntries(Object.entries(inputs).map(([name, input]) => [name, input.writes]));

    apply(['a']);
    apply(['a']);
    expectDeepEqual(writes(), {a: [false], b: [], c: []});
    apply([]);
    expectDeepEqual(writes(), {a: [false, true], b: [], c: []});
    apply([], {c: {hidden: true}});
    apply(['c'], {c: {hidden: true}});
    expectDeepEqual(writes(), {a: [false, true], b: [], c: [false]});
    applyInputsVisibility(Object.keys(inputs), (name) => inputs[name], {});
    expectDeepEqual(writes(), {a: [false, true], b: [], c: [false, true]});
  });

  test('Input form hides toggled inputs, driver meta wins', async () => {
    await awaitWebComponents();
    const root = ui.div();
    document.body.append(root);
    let form: DG.InputForm | undefined;
    const app = Vue.createApp(InputForm, {
      funcCall: DG.Func.byName('Compute2:SimpleInputs2').prepare(),
      callMeta: {c: new BehaviorSubject({hidden: true})},
      hiddenInputs: new Set(['a']),
      onFormReplaced: (f: DG.InputForm | undefined) => form = f,
    });
    app.mount(root);
    try {
      await awaitCheck(() => !!form?.getInput('c') && !form.getInput('c').visible, 'form not rendered', 10000);
      expect(form!.getInput('a').visible, false);
      expect(form!.getInput('b').visible, true);
    } finally {
      app.unmount();
      root.remove();
    }
  });

  test('Validation icon shows a link for workflow-set values', async () => {
    await awaitWebComponents();
    const icon = document.createElement('dg-validation-icon') as any;
    document.body.append(icon);
    const iconName = () => icon.querySelector('i')?.className.match(/fa-([a-z-]+)/)?.[1];
    try {
      for (const [status, expected] of [
        [consistency('restricted'), 'link'],
        [consistency('disabled'), 'link'],
        [consistency('restricted', true), 'handshake-slash'],
        [consistency('info'), undefined],
      ] as const) {
        icon.validationStatus = {consistency: status};
        await awaitCheck(() => iconName() === expected,
          `${status.restriction}: expected ${expected}, got ${iconName()}`, 3000);
      }
    } finally {
      icon.remove();
    }
  });
});
