import {Form, TextInput, button} from '@datagrok-libraries/u2';

export function newTemplateForm(onSave: (values: {name: string, key: string}) => unknown): Form {
  const form = new Form();
  form.build(() => {
    const name = new TextInput({label: 'Name', required: true});
    const key = new TextInput({label: 'Key', required: true});
    key.addValidator((k) => k === '' || /^[A-Z]{2,4}$/.test(k) ? null : '2 to 4 capital letters, like KHL');
    form.addAll([name, key]).addButtons((row) => row.append(button('SAVE',
      () => onSave({name: name.peek(), key: key.peek()}), {primary: true, disabled: form.problem})));
  });
  return form;
}
