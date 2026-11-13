/* What an app built over a platform TableView needs besides its table: rows edited in place, an ordered pick,
   confirmDialog, the dialog options and the required asterisk — and, on the platform, the TableView's chrome: a
   ribbon group owned by the view, its signals, and a dialog pinned to it. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {AccessGrid, BoolInput, ChecklistInput, ChoiceInput, Dialog, Form, RecordsInput, Scope, Signal, TextInput,
  button, computed, confirmDialog, divH, divV, h3, signal, span} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';
import {readout} from './common';

interface Field {
  name: string;
  type: string;
  required: boolean;
}

/** Fields (Name / Type / Required) and orderable stages, as a template editor has them, and the functions to
 * run in order. */
export function recordsPage(): HTMLElement {
  const fields = new RecordsInput<Field>({
    label: 'Fields', name: 'fields', addText: 'Add field',
    columns: [
      {name: 'name', header: 'Name', input: (v) => new TextInput({inline: true, value: v ?? '', required: true})},
      {name: 'type', header: 'Type', width: '120px', input: (v) => new ChoiceInput({inline: true, nullable: false,
        items: ['String', 'Number', 'Bool', 'Date', 'Molecule'], value: v ?? 'String'})},
      {name: 'required', header: 'Required', width: '64px',
        input: (v) => new BoolInput({inline: true, value: v === true})},
    ],
    newItem: () => ({name: '', type: 'String', required: false}),
    value: [{name: 'Target', type: 'String', required: true}, {name: 'Dose', type: 'Number', required: false}],
    rowName: (_, i) => `Field ${i + 1}`,
  });
  fields.addValidator((rows) => new Set(rows.map((r) => r.name)).size < rows.length ? 'Names must differ' : null);
  const stages = new RecordsInput<{name: string}>({
    label: 'Stages', name: 'stages', addText: 'Add stage', orderable: true, required: true,
    columns: [{name: 'name', header: 'Name',
      input: (v) => new TextInput({inline: true, value: v ?? '', required: true})}],
    value: [{name: 'Design'}, {name: 'Synthesis'}, {name: 'Assay'}],
    rowName: (_, i) => `Stage ${i + 1}`,
  });
  const run = new ChecklistInput({label: 'Run', name: 'run', orderable: true,
    items: ['Descriptors', 'Toxicity risks', 'Solubility', 'Structural alerts'],
    value: ['Toxicity risks', 'Descriptors']});
  const form = new Form().addAll([fields, stages, run]);
  return divV([
    span('Rows edited in place (RecordsInput: a column per field, add, remove, drag or Alt+↑/↓ to reorder) and a ' +
      'pick whose order matters (ChecklistInput). A failing cell names its row in the verdict.', 'u2demo-hint'),
    form,
    readout('stages', computed(() => stages.value.map((s) => s.name).join(', '))),
    readout('run', computed(() => run.value.join(', '))),
    readout('order', computed(() => run.order.value.join(', '))),
    readout('problem', computed(() => form.validity.value ?? 'none')),
  ], 'u2demo-page');
}

/** confirmDialog() with work that takes a moment, a sized and a resizable dialog, the first focus past a read-only
 * field, required asterisks and the access grid's headers. */
export function dialogsPage(status: Signal<string>): HTMLElement {
  const scope = Scope.ambient!;
  const deleted = signal(0);
  const remove = button('Delete campaign...', () => void confirmDialog({
    title: 'Delete campaign', name: 'delete-campaign', confirmText: 'Delete', destructive: true,
    message: 'Are you sure you want to delete campaign BTK covalent follow-up (KHL-3)?',
    detail: 'Its 260 molecules, layout and history are removed for everyone. This cannot be undone.',
    onConfirm: () => new Promise((resolve) => scope.setTimeout(() => resolve(undefined), 300)),
  }).then((yes) => {
    status.value = yes ? 'Campaign KHL-3 deleted' : 'Delete cancelled';
    if (yes)
      deleted.value++;
  }));
  const sized = button('Sized dialog', () => {
    const body = divV([span('A fixed height: the body fills it and the footer stays at the bottom.')]);
    new Dialog('Progress tracker', {name: 'sized-dialog'}).add(body).onOK(() => {}).show({width: 520, height: 420});
  });
  const resizable = button('Resizable dialog', () => {
    new Dialog('Calculate additional properties', {name: 'resizable-dialog'})
      .add(span('Drag the bottom-right corner: the body takes the size.')).onOK(() => {})
      .show({width: 480, height: 320, resizable: true});
  });
  const submit = button('Submit...', () => {
    const name = new TextInput({label: 'Name', value: 'BTK covalent follow-up', readOnly: true});
    const code = new TextInput({label: 'Code', value: 'KHL-3', readOnly: true});
    const state = new TextInput({label: 'Status', value: 'Submitted'});
    new Dialog('Submit', {name: 'submit-dialog'}).add(new Form().addAll([name, code, state]))
      .onOK(() => status.value = `Submitted as ${state.peek()}`, {text: 'Submit'}).show();
  });
  const required = new Form().addAll([
    new TextInput({label: 'Name', name: 'template-name', required: true}),
    new TextInput({label: 'Key', name: 'template-key', required: true}),
    new TextInput({label: 'Description', name: 'template-description'}),
  ]);
  const access = new AccessGrid({label: 'Access', name: 'access', capabilities: [{name: 'view', label: 'Can View'},
    {name: 'edit', label: 'Can Edit'}], principals: ['All users', 'Chemists'],
  value: [{principal: 'Chemists', can: {view: true, edit: true}}]});
  return divV([
    span('confirmDialog() asks before what cannot be undone: CANCEL and the action as a verb. Dialogs with a fixed ' +
      'or a user-set height let the body fill them; the first focus skips read-only fields.', 'u2demo-hint'),
    divH([remove, sized, resizable, submit], 'u2demo-row'),
    readout('deleted', computed(() => String(deleted.value))),
    h3('Required fields'),
    required,
    h3('Access'),
    access,
  ], 'u2demo-page');
}

/** Opens a TableView over a small campaign with a u2 ribbon group owned by the view, the view's signals as a
 * readout here and in the view's status bar, and a progress tracker pinned to the view. */
export function tableChromePage(status: Signal<string>): HTMLElement {
  const opened = signal('none');
  const row = signal('');
  const open = button('Open the campaign view', async () => {
    const df = DG.DataFrame.fromCsv('Molecule,Stage,MolWt\nCCO,Design,46\nc1ccccc1,Synthesis,78\nCC(=O)O,Assay,60');
    df.name = 'KHL-3';
    const home = grok.shell.v;
    // as an app opens a campaign: after its table is read — and after the click, whose end makes the view under
    // the pointer current again
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    // the session lives as long as the view, not as this page
    const host = ViewHost.of(tv);
    const signals = host.signals;
    const changes = signal<string | null>(null);
    const canEdit = signal(true);
    const done = signal('');
    const say = (action: string) => {
      done.value = action;
      status.value = action;
    };
    const line = computed(() => `Row ${signals.currentRowIdx.value + 1} of ${signals.rowCount.value}` +
      (done.value === '' ? '' : ` · ${done.value}`));
    host.effect(() => {
      row.value = line.value;
    });
    host.own(() => opened.value = 'closed');
    host.build(() => host.ribbon({
      hideSave: true, name: 'demo-campaign-actions', status: line,
      breadcrumbs: {parent: home, onClick: () => tv.close()},
      actions: [
        {icon: 'sync', tooltip: 'Refresh Campaign', run: () => {
          changes.value = null;
          say('Refreshed');
        }, hint: changes, hintDetail: 'New changes'},
        {icon: 'plus', tooltip: 'Add new row', run: () => void df.rows.addNew(null, true)},
        {icon: 'code-merge', tooltip: 'Merge another table into the campaign',
          text: 'Merge another table into the campaign...', run: () => say('Merge')},
        {icon: 'wrench', tooltip: 'Calculate additional properties', text: 'Calculate additional properties...',
          run: () => say('Calculate')},
        {icon: 'users', tooltip: 'Edit campaign permissions', text: 'Edit campaign permissions...',
          run: () => say('Permissions'), visible: canEdit},
        {icon: 'th-large', tooltip: 'Apply template layout (Static)', run: () => say('Layout')},
        {icon: 'bell', tooltip: 'Simulate new changes', run: () => {
          changes.value = 'New changes available. Click to refresh campaign.';
        }},
      ],
      buttons: [
        {text: 'Progress tracker', run: () => void host.show(new Dialog('Progress tracker', {name: 'demo-tracker'})
          .add(span('Hidden while another view is current; closed with the campaign view.')).onOK(() => {}),
        {modal: true})},
        {text: 'Submit...', run: () => say('Submit'), cls: 'hit-design-submit-button'},
      ],
    }));
    opened.value = 'open';
  });
  return divV([
    span('A TableView owned by u2: the package group sits after a rule at the ribbon\'s end, pinned while the ribbon ' +
      'overflows, its icons folding into "…" when there is no room; the crumb goes back here; the project Save is ' +
      'left off; the status bar and this readout follow the view\'s current row.', 'u2demo-hint'),
    divH([open], 'u2demo-row'),
    readout('campaign view', opened),
    readout('current', row),
  ], 'u2demo-page');
}
