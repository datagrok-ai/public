import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Dialog, Form, TextInput, computed, iconButton, signal} from '@datagrok-libraries/u2';

export function mergeDialog(onMerge: (file: DG.FileInfo) => unknown): Dialog {
  const dialog = new Dialog('Merge with table');
  const chosen = signal<DG.FileInfo | null>(null);
  dialog.add(() => {
    const path = new TextInput({label: 'Source file', placeholder: 'System:DemoFiles/file.csv',
      onInput: async (text) => {
        chosen.value = null;
        if (await grok.dapi.files.exists(text)) {
          const files = await grok.dapi.files.list(text.substring(0, text.lastIndexOf('/') + 1));
          chosen.value = files.find((f) => f.fullPath === text) ?? null;
        }
      }});
    path.addOptions(iconButton('folder-open', () => {
      const picker = document.createElement('input');
      picker.type = 'file';
      picker.onchange = async () => {
        const file = picker.files?.[0];
        if (!file)
          return;
        chosen.value = DG.FileInfo.fromBytes(file.name, new Uint8Array(await file.arrayBuffer()));
        path.value = file.name;
      };
      picker.click();
    }, {tooltip: 'Open file'}));
    return new Form().add(path);
  });
  dialog.onOK(() => onMerge(chosen.peek()!), {text: 'MERGE'})
    .okDisabled(computed(() => chosen.value === null ? 'Choose a file' : null));
  return dialog.show();
}
