import {Form, TextArea, TextInput, button} from '@datagrok-libraries/u2';

export function templateEditor(onCreate: (values: {name: string, key: string, description: string}) => unknown,
  onCancel: () => void): Form {
  const form = new Form();
  form.root.style.maxWidth = '640px';
  form.build(() => {
    const name = new TextInput({label: 'Name'});
    const key = new TextInput({label: 'Key'});
    const keyEditor = key.root.querySelector<HTMLElement>('input')!;
    keyEditor.style.width = '120px';
    keyEditor.style.flexGrow = '0';
    const description = new TextArea({label: 'Description'});
    form.addAll([name, key, description]).addButtons((row) => {
      row.style.gap = '8px';
      row.append(
        button('CANCEL', () => onCancel()),
        button('CREATE', () => onCreate({name: name.peek(), key: key.peek(),
          description: description.peek()}), {primary: true}));
    });
  });
  return form;
}
