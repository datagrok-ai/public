import {BoolInput, ChoiceInput, Control, Form, RecordsInput, TextInput, button, component, notify}
  from '@datagrok-libraries/u2';

interface Field {
  name: string;
  type: string;
  required: boolean;
}

const TYPES = ['Text', 'Number', 'Date', 'Molecule'];

export function templateFields(initial: Field[], save: (fields: Field[]) => void): Control {
  return component(() => {
    const fields = new RecordsInput<Field>({
      label: 'Fields', orderable: true, addText: 'Add field', required: true,
      columns: [
        {name: 'name', header: 'Name', input: (v) => new TextInput({inline: true, value: v ?? '', required: true})},
        {name: 'type', header: 'Type', input: (v) => new ChoiceInput({inline: true, items: TYPES, value: v ?? 'Text'})},
        {name: 'required', header: 'Required', input: (v) => new BoolInput({inline: true, value: v ?? false})},
      ],
      newItem: () => ({name: '', type: 'Text', required: false}),
      value: initial,
      rowName: (field, i) => field.name || `Field ${i + 1}`,
    });
    fields.addValidator((rows) => new Set(rows.map((r) => r.name)).size < rows.length ? 'Field names must differ' : null);
    const form = new Form().add(fields);
    return form.addButtons((row) => row.append(
      button('Save template', () => save(fields.value), {primary: true, disabled: form.problem})));
  }, {name: 'template-fields'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = templateFields([{name: 'Compound', type: 'Molecule', required: true},
    {name: 'IC50', type: 'Number', required: false}], (fields) => notify.info(`Saved ${fields.length} fields`));
  host.append(demo.root);
  return demo;
}
