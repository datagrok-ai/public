import {BoolInput, ChoiceInput, Control, DateInput, Form, Format, NumberInput, SliderInput, TextArea, TextInput, component,
  computed, statusLine} from '@datagrok-libraries/u2';

export function dosingForm(): Control {
  return component(() => {
    const name = new TextInput({label: 'Name', value: 'Aspirin'});
    const dose = new NumberInput({label: 'Dose', value: 50, min: 0, max: 1000, postfix: 'mg'});
    const route = new ChoiceInput({label: 'Route', items: ['Oral', 'IV', 'Topical'], value: 'Oral'});
    const start = new DateInput({label: 'Start', value: new Date(2026, 9, 1)});
    const confidence = new SliderInput({label: 'Confidence', min: 0, max: 100, value: 80});
    const active = new BoolInput({label: 'Active', value: true, switch: true});
    const notes = new TextArea({label: 'Notes', placeholder: 'Anything the reviewer should know',
      disabled: computed(() => active.value ? null : 'Switch Active on to add notes')});
    const form = new Form().addAll([name, dose, route, start, confidence, active, notes]);
    const summary = computed(() => `${name.value}: ${dose.value} mg ${route.value?.toLowerCase()}, ` +
      `from ${Format.value(start.value)}, ${confidence.value}% sure`);
    return form.addElement(statusLine(form.shownProblem, summary));
  }, {name: 'dosing-form'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = dosingForm();
  host.append(demo.root);
  return demo;
}
