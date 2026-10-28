import {Control, component, registerAll, renderSpec} from '@datagrok-libraries/u2';
import type {Spec} from '@datagrok-libraries/u2';

const DOSE: Spec = {
  $schema: 'dg-ui/1',
  root: {tag: 'u2-div-v', children: [
    {tag: 'h3', bind: {text: '$.compound'}},
    {tag: 'u2-form', children: [
      {tag: 'u2-text-input', name: 'compound', props: {label: 'Compound', value: 'Aspirin'}},
      {tag: 'u2-radio-input', name: 'unit', props: {label: 'Units', items: ['mg', 'mg/kg'], value: 'mg', buttons: true}},
      {tag: 'u2-number-input', props: {label: 'Dose', value: 100, min: 0}, bind: {postfix: '$.unit'}},
      {tag: 'u2-bool-input', name: 'advanced', props: {label: 'More options', switch: true}},
      {tag: 'u2-text-area', props: {label: 'Notes'}, bind: {visible: '$.advanced'}},
    ]},
  ]},
};

export function doseSpec(): Control {
  registerAll();
  return component(() => renderSpec(DOSE), {name: 'dose-spec'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = doseSpec();
  host.append(demo.root);
  return demo;
}
