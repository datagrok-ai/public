import {Form, TextArea, TextInput, button} from '@datagrok-libraries/u2';
import './good.css';

export function templateEditor(onCreate: (values: {name: string, key: string, description: string}) => unknown,
  onCancel: () => void): Form {
  const form = new Form();
  form.root.classList.add('template-editor');
  form.build(() => {
    const name = new TextInput({label: 'Name', fill: true});
    const key = new TextInput({label: 'Key', width: 120});
    const description = new TextArea({label: 'Description', fill: true});
    form.addAll([name, key, description]).addButtons((row) => row.append(
      button('CANCEL', () => onCancel()),
      button('CREATE', () => onCreate({name: name.peek(), key: key.peek(),
        description: description.peek()}), {primary: true})));
  });
  return form;
}
