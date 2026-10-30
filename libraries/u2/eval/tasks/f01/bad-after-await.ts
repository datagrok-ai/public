import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChoiceInput, Dialog, Form, TextInput, computed, effect, iconButton, signal, untracked}
  from '@datagrok-libraries/u2';
import {FileInput, Platform, fileInput, ViewHost} from '@datagrok-libraries/u2/dg';

type Campaign = {name: string, filePath: string, status: string};

interface CampaignService {
  rename(name: string): Promise<void>;
  download(): void;
  submit(): Promise<void>;
  readTable(file: DG.FileInfo): Promise<DG.DataFrame>;
  merge(file: DG.FileInfo, table: DG.DataFrame, moleculeColumn: string): Promise<void>;
  refresh(): Promise<void>;
}

function submitDialog(campaign: Campaign, service: CampaignService): Dialog {
  const dialog = new Dialog('Submit campaign');
  dialog.add(() => {
    const editing = signal(false);
    const finish = (save: boolean) => {
      if (save) {
        campaign.name = name.peek();
        service.rename(campaign.name).catch((e) => grok.shell.error(`Rename failed: ${e}`));
      }
      else
        name.reset(campaign.name);
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
    path.addOptions(iconButton('download', () => service.download(), {tooltip: 'Download'}));
    const status = new TextInput({label: 'Status', value: campaign.status, readOnly: true});
    return new Form().addAll([name, path, status]);
  });
  dialog.onOK(() => service.submit(), {text: 'SUBMIT'});
  return dialog;
}

const isMolecule = (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE;

function mergeDialog(service: CampaignService): Dialog {
  const dialog = new Dialog('Merge with table');
  const table = signal<DG.DataFrame | null>(null);
  const molecules = signal<ChoiceInput | null>(null);
  let form!: Form;
  let source!: FileInput;
  const read = async (file: DG.FileInfo | null) => {
    molecules.peek()?.root.remove();
    molecules.value = null;
    table.value = null;
    if (file === null)
      return;
    const df = await service.readTable(file);
    const names = Platform.listedColumns(df, isMolecule).map((c) => c.name);
    const input = new ChoiceInput({label: 'Molecules', items: names, value: names[0] ?? null});
    form.add(input);
    table.value = df;
    molecules.value = input;
  };
  dialog.add(() => {
    source = fileInput('Source file', {required: true, onChanged: (file) => void read(file)});
    return form = new Form().add(source);
  });
  dialog.onOK(() => service.merge(source.peek()!, table.peek()!, molecules.peek()!.peek()!),
    {text: 'MERGE'})
    .okDisabled(computed(() => source.value === null ? 'Choose the file to merge' :
      molecules.value === null ? 'Reading the file…' : null));
  return dialog;
}

export function campaignActions(tv: DG.TableView, campaign: Campaign, service: CampaignService): void {
  const record = {...campaign};
  ViewHost.of(tv).ribbon({
    name: 'campaign-actions',
    actions: [
      {icon: 'plus', tooltip: 'Add new row', run: () => {
        tv.dataFrame.rows.addNew(null, true);
      }},
      {icon: 'sync', tooltip: 'Refresh Campaign', run: () => service.refresh()},
      {icon: 'code-merge', tooltip: 'Merge with table', text: 'Merge with table...', run: () => {
        ViewHost.of(tv).show(mergeDialog(service));
      }},
    ],
    buttons: [{text: 'Submit...', run: () => {
      ViewHost.of(tv).show(submitDialog(record, service));
    }}],
  });
}
