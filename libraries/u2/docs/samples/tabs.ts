import {ChoiceInput, Control, Form, NumberInput, TabStrip, TextArea, TextInput, Toolbar, component, notify}
  from '@datagrok-libraries/u2';

export function experimentEditor(): Control {
  return component(() => {
    const general = new Form().addAll([new TextInput({label: 'Name', required: true}),
      new ChoiceInput({label: 'Owner', items: ['Chemistry', 'Biology', 'DMPK'], value: 'Chemistry'})]);
    const assay = new Form().addAll([new NumberInput({label: 'Concentration', min: 0, required: true, postfix: 'µM'}),
      new NumberInput({label: 'Replicates', mode: 'int', min: 1, max: 12, value: 3})]);
    const notes = new Form().add(new TextArea({label: 'Notes', rows: 4}));
    const tabs = new TabStrip({tabs: [
      {id: 'general', title: 'General', content: general},
      {id: 'assay', title: 'Assay', content: assay},
      {id: 'notes', title: 'Notes', content: notes},
    ]});
    const save = () => {
      if ([general, assay, notes].map((form) => form.validate()).every((valid) => valid))
        notify.info('Saved');
    };
    return [tabs, new Toolbar({ariaLabel: 'Experiment'}).addButton('Save', save, {icon: 'save'})];
  }, {name: 'experiment-editor'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = experimentEditor();
  host.append(demo.root);
  return demo;
}
