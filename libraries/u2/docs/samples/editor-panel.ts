import {ChoiceInput, ColorInput, Control, Form, Section, TextArea, TextInput, button, component, computed, notify,
  statusLine} from '@datagrok-libraries/u2';

export function monomerEditor(): Control {
  return component(() => {
    const symbol = new TextInput({label: 'Symbol', value: 'Aib', required: true});
    const name = new TextInput({label: 'Name', value: 'α-Aminoisobutyric acid'});
    const polymer = new ChoiceInput({label: 'Polymer', items: ['PEPTIDE', 'RNA', 'CHEM'], value: 'PEPTIDE'});
    const type = new ChoiceInput({label: 'Type', items: ['Backbone', 'Branch', 'Terminal'], value: 'Backbone'});
    const colors = [new ColorInput({label: 'Line', value: '#2083D5'}),
      new ColorInput({label: 'Background', value: '#ECEFF2'}), new ColorInput({label: 'Text', value: '#4D5261'})];
    const notes = new TextArea({label: 'Notes', placeholder: 'Where it comes from, how to use it'});
    const details = new Section({title: 'Details', hint: computed(() => notes.value || 'Colors and notes'),
      validity: notes.validity, expanded: false, chevron: 'end', persistKey: 'build-with-u2.monomer.details'});
    const form = new Form()
      .addColumns([symbol, name, polymer, type], {minWidth: 176})
      .addElement(details.root)
      .addGroup('Colors', colors, {captionsAfter: true, host: details.body})
      .add(notes, details.body);
    return form.addButtons((row) => row.append(
      button('Save', () => notify.info(`Saved ${symbol.value}`), {primary: true, disabled: form.problem}),
      statusLine(form.shownProblem, computed(() => `Saves ${symbol.value} into the ${polymer.value} library`)),
    ));
  }, {name: 'monomer-editor'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = monomerEditor();
  host.append(demo.root);
  return demo;
}
