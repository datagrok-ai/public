import {Accordion, Control, Form, NumberInput, TextArea, component, computed, tableFromMap} from '@datagrok-libraries/u2';

export function compoundPanel(): Control {
  return component(() => {
    const dose = new NumberInput({label: 'Dose', value: 10, min: 0, max: 100, postfix: 'mg/kg'});
    const notes = new TextArea({label: 'Notes', placeholder: 'At most 80 characters'});
    notes.addValidator((v) => v.length > 80 ? `${v.length - 80} characters too long` : null);
    const accordion = new Accordion();
    accordion.addPane('Compound', () => tableFromMap({'Name': 'Aspirin', 'Formula': 'C9H8O4', 'MW': '180.16 g/mol'}),
      {expanded: true});
    accordion.addPane('Dosing', new Form().add(dose).root,
      {hint: computed(() => `${dose.value ?? '?'} mg/kg`), persistKey: 'build-with-u2.dosing'});
    accordion.addPane('Notes', new Form().add(notes).root,
      {hint: computed(() => notes.value ? 'Has notes' : 'None'), validity: notes.validity});
    return accordion;
  }, {name: 'compound-panel'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = compoundPanel();
  host.append(demo.root);
  return demo;
}
