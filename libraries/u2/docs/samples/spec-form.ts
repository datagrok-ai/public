import {Control, SpecContext, component, computed, notify, registerAll, renderSpec, signal} from '@datagrok-libraries/u2';
import type {Spec} from '@datagrok-libraries/u2';

const DOSING: Spec = {
  $schema: 'dg-ui/1',
  root: {tag: 'u2-form', children: [
    {tag: 'u2-text-input', props: {label: 'Compound', required: true}, bind: {value: '$.compound'}},
    {tag: 'u2-number-input', props: {label: 'Dose', min: 0, postfix: 'mg'}, bind: {value: '$.dose'}},
    {tag: 'u2-choice-input', props: {label: 'Route', items: ['Oral', 'IV', 'Topical']}, bind: {value: '$.route'}},
    {tag: 'u2-form-buttons', children: [{tag: 'u2-button', props: {text: 'Save', primary: true},
      bind: {disabled: '$.saveProblem'}, on: {click: 'cmd:save'}}]},
  ]},
};

export function dosingForm(): Control {
  registerAll();
  const compound = signal('Aspirin');
  const dose = signal<number | null>(100);
  const route = signal<string | null>('Oral');
  const saveProblem = computed(() => compound.value.trim() ? null : 'Enter the compound');
  return component(() => renderSpec(DOSING, new SpecContext({
    data: {compound, dose, route, saveProblem},
    commands: {save: () => notify.info(`Saved ${compound.value}: ${dose.value} mg ${route.value?.toLowerCase()}`)},
  })), {name: 'dosing-form'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = dosingForm();
  host.append(demo.root);
  return demo;
}
