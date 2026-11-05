import {Form, TextInput, button, computed, effect} from '@datagrok-libraries/u2';

export function newTemplateForm(onSave: (values: {name: string, key: string}) => unknown): Form {
  const form = new Form();
  form.build(() => {
    const name = new TextInput({label: 'Name', required: true});
    const key = new TextInput({label: 'Key', required: true});
    key.addValidator((k) => k === '' || /^[A-Z]{2,4}$/.test(k) ? null : '2 to 4 capital letters, like KHL');
    const invalid = computed(() => form.validity.value !== null);
    const save = button('SAVE', () => onSave({name: name.peek(), key: key.peek()}), {primary: true});
    effect(() => save.disabled = invalid.value);
    form.addAll([name, key]).addButtons((row) => row.append(save));
  });
  return form;
}
