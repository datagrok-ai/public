import {Control, Form, Page, TextInput, button, component, computed, iconButton, notify, statusLine}
  from '@datagrok-libraries/u2';

export function templateEditor(): Control {
  return component(() => {
    const name = new TextInput({label: 'Name', required: true});
    const key = new TextInput({label: 'Key', required: true});
    const description = new TextInput({label: 'Description'});
    const details = new Form().addAll([name, key]);
    const more = new Form().add(description);
    const create = button('Create', () => notify.info(`Created ${name.peek()}`),
      {primary: true, disabled: details.problem});
    return new Page({maxWidth: 820,
      sections: [
        {title: 'Details', content: details},
        {title: 'Description', content: more,
          actions: [iconButton('eraser', () => description.value = '', {tooltip: 'Clear the description'})]},
      ],
      footer: [statusLine(details.shownProblem, computed(() => `Creates the template ${key.value || '…'}`)), create]});
  }, {name: 'template-editor'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const editor = templateEditor();
  host.append(editor.root);
  return editor;
}
