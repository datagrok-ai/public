import {ChecklistInput, Control, Form, component, computed, statusLine} from '@datagrok-libraries/u2';

const STEPS = ['Calculate descriptors', 'Filter by Lipinski', 'Predict solubility', 'Dock against the target', 'Score and rank'];

export function pipelineSteps(chosen: string[]): Control {
  return component(() => {
    const steps = new ChecklistInput({label: 'Pipeline', items: STEPS, orderable: true, required: true, value: chosen});
    const plan = computed(() => steps.value.map((s, i) => `${i + 1}. ${s}`).join(' → '));
    const form = new Form().add(steps);
    return form.addElement(statusLine(form.shownProblem, plan));
  }, {name: 'pipeline-steps'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = pipelineSteps(['Calculate descriptors', 'Score and rank']);
  host.append(demo.root);
  return demo;
}
