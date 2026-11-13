import {BoolInput, Control, Form, MasterDetail, NumberInput, Signal, component, computed, div, signal}
  from '@datagrok-libraries/u2';
import './item-settings.css';

interface Calculation {
  key: string;
  name: string;
  on: Signal<boolean>;
}

export function calculations(items: Calculation[]): Control {
  return component(() => div([new MasterDetail({items, keyOf: (c) => c.key, name: 'calculations',
    row: (c) => ({label: c.name, check: c.on, hint: computed(() => c.on.value ? 'runs' : '')}),
    detail: (c) => new Form().addAll([new NumberInput({label: 'Threshold', name: `${c.key}-threshold`, value: 5}),
      new BoolInput({label: 'Run again on open', name: `${c.key}-rerun`})]),
    disabled: (c) => computed(() => c.on.value ? null : `Check ${c.name} to change its settings`)})],
  'calculations'), {name: 'calculations-editor'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const editor = calculations(['Solubility', 'Toxicity', 'Lipophilicity'].map((name, i) =>
    ({key: name.toLowerCase(), name, on: signal(i === 0)})));
  host.append(editor.root);
  return editor;
}
