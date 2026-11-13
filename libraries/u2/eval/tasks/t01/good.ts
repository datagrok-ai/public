import {Dialog, Form, TextInput, computed, effect, iconButton, signal} from '@datagrok-libraries/u2';

export function submitDialog(campaign: {name: string, filePath: string, status: string},
  actions: {rename(name: string): unknown, download(): void, submit(): unknown}): Dialog {
  const dialog = new Dialog('Submit campaign');
  dialog.add(() => {
    const editing = signal(false);
    let saved = campaign.name;
    const finish = (save: boolean) => {
      if (save) {
        saved = name.peek();
        actions.rename(saved);
      }
      else
        name.reset(saved);
      editing.value = false;
    };
    const name: TextInput = new TextInput({label: 'Name', value: campaign.name,
      readOnly: computed(() => !editing.value),
      onEnter: () => editing.peek() && finish(true),
      onEscape: () => editing.peek() ? finish(false) : dialog.cancelByKey()});
    const edit = iconButton('pencil', () => {
      editing.value = true;
      name.focus();
    }, {tooltip: 'Edit name'});
    const save = iconButton('check', () => finish(true), {tooltip: 'Save name'});
    const cancel = iconButton('times', () => finish(false), {tooltip: 'Cancel'});
    effect(() => {
      edit.hidden = editing.value;
      save.hidden = cancel.hidden = !editing.value;
    });
    for (const el of [edit, save, cancel])
      name.addOptions(el);
    const path = new TextInput({label: 'File Path', value: campaign.filePath, readOnly: true});
    path.addOptions(iconButton('download', () => actions.download(), {tooltip: 'Download'}));
    const status = new TextInput({label: 'Status', value: campaign.status, readOnly: true});
    return new Form().addAll([name, path, status]);
  });
  dialog.onOK(() => actions.submit(), {text: 'SUBMIT'});
  return dialog.show();
}
