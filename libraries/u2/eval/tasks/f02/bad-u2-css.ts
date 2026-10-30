import * as DG from 'datagrok-api/dg';
import {BoolInput, Control, Form, Input, NumberInput, button, computed, span} from '@datagrok-libraries/u2';
import {tableInput} from '@datagrok-libraries/u2/dg';
import './bad-u2-css.css';

type Settings = {[func: string]: {[param: string]: unknown}};

interface ComputeStore {
  load(): Promise<{functions: DG.Func[], saved: Settings}>;
  save(settings: Settings): Promise<void>;
}

type Entry = {func: DG.Func, check: BoolInput, fields: [string, Input<any>][]};

function settingsForm(functions: DG.Func[], saved: Settings, store: ComputeStore): Form {
  const form = new Form();
  const entries: Entry[] = [];
  for (const func of functions) {
    const caption = func.friendlyName ?? func.name;
    const values = saved[func.name] ?? {};
    const check = new BoolInput({label: caption, value: func.name in saved});
    const on = computed(() => check.value);
    const table = func.inputs.find((p) => p.propertyType === DG.TYPE.DATA_FRAME);
    const fields = func.inputs.filter((p) => p !== table).map((p): [string, Input<any>] => [p.name,
      p.propertyType === DG.TYPE.DATA_FRAME ?
        tableInput(p.name, {value: (values[p.name] as string | undefined) ?? null, disabled: computed(() => on.value ? null : `Check ${caption} to change its settings`)}) :
        new NumberInput({label: p.name, value: (values[p.name] as number | undefined) ?? p.defaultValue ?? null,
          disabled: computed(() => on.value ? null : `Check ${caption} to change its settings`)})]);
    form.add(check).addAll(fields.map(([, input]) => input));
    if (fields.length > 0)
      form.addElement(span(computed(() => on.value ? '' : `Check ${caption} to change its settings`)));
    entries.push({func, check, fields});
  }
  const settings = (): Settings => Object.fromEntries(entries.filter((e) => e.check.peek())
    .map((e) => [e.func.name, Object.fromEntries(e.fields.map(([name, input]) => [name, input.peek()]))]));
  const none = computed(() => entries.every((e) => !e.check.value));
  form.addButtons((row) => row.append(button('SAVE', () => store.save(settings()), {primary: true,
    disabled: computed(() => none.value ? 'Check the functions to run first' : null)})));
  return form;
}

export function computeSettings(table: DG.DataFrame, store: ComputeStore): Control {
  const pane = new Control();
  pane.root.classList.add('compute-settings');
  pane.root.append(span('Loading…'));
  store.load().then(({functions, saved}) => {
    if (pane.scope.isDisposed)
      return;
    pane.build(() => pane.root.replaceChildren(settingsForm(functions, saved, store).root));
  }, (e) => pane.root.replaceChildren(span(`Could not load the compute settings: ${e}`)));
  return pane;
}
