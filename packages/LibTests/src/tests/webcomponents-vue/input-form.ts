import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BehaviorSubject} from 'rxjs';
import {category, test, expect, waitUntil, delay} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {InputForm} from '@datagrok-libraries/webcomponents-vue';
import {
  applyInputsVisibility, injectInputBaseStatus, isInputEnabled,
} from '@datagrok-libraries/webcomponents-vue/src/InputForm/utils';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {awaitWebComponents} from '../../test-utils';
import {ADVICE, consistency, findLink, openPopover, scriptWith} from '../webcomponents/fixtures';

async function mountInputForm(props: () => Record<string, any>) {
  await awaitWebComponents();
  const root = ui.div();
  document.body.append(root);
  let form: DG.InputForm | undefined;
  const app = Vue.createApp({
    render: () => Vue.h(InputForm as Vue.Component,
      {...props(), onFormReplaced: (f: DG.InputForm | undefined) => form = f}),
  });
  app.mount(root);
  return {
    input: (name: string) => form!.getInput(name),
    awaitForm: (name: string) => waitUntil(() => !!form?.getInput(name), 'form not rendered', 10000),
    dispose: () => {
      app.unmount();
      root.remove();
    },
  };
}

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

const marks = (input: DG.InputBase) => [input.input, input.captionLabel].map((el) =>
  ['d4-invalid', 'd4-partially-invalid'].filter((cls) => el.classList.contains(cls)).join(' '));

category('WebComponentsVue: InputForm', () => {
  test('mutateOnly locks restricted and info dataframe inputs only', async () => {
    const [df, plain, x] = scriptWith('dataframe df {mutateOnly: true}', 'dataframe plain',
      'double x {mutateOnly: true}').inputs;
    const enabled = (prop: DG.Property, isReadonly = false) => Object.fromEntries(
      (['disabled', 'restricted', 'info', 'none', undefined] as const).map((restriction) =>
        [restriction ?? 'unlinked', isInputEnabled(prop, restriction && consistency(restriction), isReadonly)]));
    expectDeepEqual(enabled(df), {disabled: false, restricted: false, info: false, none: true, unlinked: true});
    expectDeepEqual(enabled(plain), {disabled: false, restricted: true, info: true, none: true, unlinked: true});
    expectDeepEqual(enabled(x), {disabled: false, restricted: true, info: true, none: true, unlinked: true});
    expectDeepEqual(enabled(plain, true),
      {disabled: false, restricted: false, info: false, none: false, unlinked: false});
  });

  test('visibility is written only when it changes', async () => {
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

  test('hidden inputs and the hidden meta both hide an input', async () => {
    const call = scriptWith('int a', 'double b', 'string c').prepare({a: 1, b: 2, c: 'x'});
    const h = await mountInputForm(() => ({
      funcCall: call,
      callMeta: {c: new BehaviorSubject({hidden: true})},
      hiddenInputs: new Set(['a']),
    }));
    try {
      await h.awaitForm('c');
      await waitUntil(() => !h.input('c').visible, 'the hidden meta did not hide c', 5000);
      expect(h.input('a').visible, false);
      expect(h.input('b').visible, true);
    } finally {
      h.dispose();
    }
  });

  test('an injected input forwards advice params with the action request', async () => {
    await awaitWebComponents();
    const input = ui.input.forProperty(DG.Property.js('x', DG.TYPE.FLOAT));
    const emitted: any[] = [];
    injectInputBaseStatus((...args: any[]) => emitted.push(args), 'x', input);
    document.body.append(input.root);
    try {
      (input as any).setStatus(ADVICE);
      findLink(await openPopover(input.root), 'Set to 10')!.click();
      expectDeepEqual(emitted, [['actionRequested', 'action-1', {val: 10}]]);
    } finally {
      input.root.remove();
    }
  });

  test('driver validation marks the editor and the caption', async () => {
    const call = scriptWith('double x').prepare({x: 1});
    const states = Vue.shallowRef<Record<string, any>>({});
    const h = await mountInputForm(() => ({funcCall: call, skipLogic: true, validationStates: states.value}));
    try {
      await h.awaitForm('x');
      await delay(500);
      states.value = {x: {errors: ['bad']}};
      await waitUntil(() => marks(h.input('x')).every((m) => m === 'd4-invalid'),
        () => `errors: ${marks(h.input('x'))}`, 3000);
      states.value = {x: {warnings: ['meh']}};
      await waitUntil(() => marks(h.input('x')).every((m) => m === 'd4-partially-invalid'),
        () => `warnings: ${marks(h.input('x'))}`, 3000);
      states.value = {};
      await waitUntil(() => marks(h.input('x')).every((m) => m === ''),
        () => `cleared: ${marks(h.input('x'))}`, 3000);
    } finally {
      h.dispose();
    }
  });

  // workflow forms (skipLogic) have no platform validation; standalone forms take no driver states
  test('driver errors present when a workflow form opens stay marked', async () => {
    const call = scriptWith('double x').prepare({x: 1});
    const states = {x: {errors: ['bad']}};
    const h = await mountInputForm(() => ({funcCall: call, skipLogic: true, validationStates: states}));
    try {
      await h.awaitForm('x');
      await delay(500);
      expectDeepEqual(marks(h.input('x')), ['d4-invalid', 'd4-invalid']);
    } finally {
      h.dispose();
    }
  });

  test('toggling isReadonly alone updates the inputs', async () => {
    const call = scriptWith('double x').prepare({x: 1});
    const readonly = Vue.ref(false);
    const states = {};
    const h = await mountInputForm(() => ({funcCall: call, isReadonly: readonly.value, validationStates: states}));
    try {
      await h.awaitForm('x');
      await waitUntil(() => h.input('x').enabled, 'the input should start enabled', 3000);
      readonly.value = true;
      await Vue.nextTick();
      await delay(100);
      expect(h.input('x').enabled, false, 'isReadonly did not disable the input');
    } finally {
      h.dispose();
    }
  });

  test('a nullable choice offers an empty option, the emptyChoice meta overrides it', async () => {
    const cases: [string, Record<string, any> | undefined, boolean][] = [
      ['string c {choices: ["a", "b"]; nullable: true}', undefined, true],
      ['string c {choices: ["a", "b"]}', undefined, false],
      ['string c {choices: ["a", "b"]; nullable: true}', {emptyChoice: false}, false],
      ['string c {choices: ["a", "b"]}', {emptyChoice: true}, true],
    ];
    for (const [annotation, meta, nullable] of cases) {
      const call = scriptWith(annotation).prepare({c: 'a'});
      const h = await mountInputForm(() => ({
        funcCall: call, callMeta: meta ? {c: new BehaviorSubject(meta)} : undefined,
      }));
      try {
        await h.awaitForm('c');
        await waitUntil(() => h.input('c').nullable === nullable,
          `${annotation} ${JSON.stringify(meta)}: expected nullable ${nullable}`, 3000);
        expect(h.input('c').value, 'a', 'the empty option must not change the value');
      } finally {
        h.dispose();
      }
    }
  });

  test('literal choices fill the items under skipLogic', async () => {
    const call = scriptWith('string c {choices: ["x", "y"]}').prepare({c: 'y'});
    const h = await mountInputForm(() => ({funcCall: call, skipInit: false, skipLogic: true}));
    try {
      await h.awaitForm('c');
      await waitUntil(() => JSON.stringify((h.input('c') as DG.ChoiceInput<string>).items) === '["x","y"]',
        () => `items: ${JSON.stringify((h.input('c') as DG.ChoiceInput<string>).items)}`, 3000);
      expect(h.input('c').value, 'y');
    } finally {
      h.dispose();
    }
  });

  test('float inputs show Float32 values without noise', async () => {
    const call = scriptWith('double x').prepare({x: Math.fround(0.2)});
    const h = await mountInputForm(() => ({funcCall: call}));
    try {
      await h.awaitForm('x');
      await waitUntil(() => (h.input('x').input as HTMLInputElement).value === '0.2',
        () => `got ${(h.input('x').input as HTMLInputElement).value}`, 3000);
    } finally {
      h.dispose();
    }
  });

  test('an exported float input matches the form, shown or not', async () => {
    const exported = async (call: DG.FuncCall) => {
      const [, wb] = await richFunctionViewReport('Excel', call.func, call, {});
      return String(wb.getWorksheet('Input scalars')!.getRow(2).getCell(2).value);
    };
    const call = scriptWith('double x').prepare({x: 0.0001});
    const beforeForm = await exported(call);
    const h = await mountInputForm(() => ({funcCall: call}));
    try {
      await h.awaitForm('x');
      const shown = (h.input('x').input as HTMLInputElement).value;
      expect(shown, '0.0001');
      expect(await exported(call), shown, 'export after the form was shown');
      expect(beforeForm, shown, 'export of a step whose form was never shown');
    } finally {
      h.dispose();
    }
  });
});
