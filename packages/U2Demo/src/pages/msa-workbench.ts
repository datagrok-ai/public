/* MSA workbench — a small self-contained app the behavioral tests in libraries/bdd drive
   (features/demo/*.feature). Automation identity comes from the u2 contract alone: every control
   stamps `data-u2`, deliberate names stamp `data-u2-name`, parts stamp `data-u2-part`; the
   feature files address the page as "run button in toolbar", "sequence column input in MSA
   dialog", "second item in aligned sequences list". */
import {
  signal,
  button, div, divH, divV, h1, h2, span,
  TextInput, ChoiceInput, NumberInput, BoolInput, Form, Toolbar, Dialog, TabStrip, VirtualList, Section, notify,
} from '@datagrok-libraries/u2';
import '../../css/u2demo-msa.css';

const SEQUENCE_COLUMNS = ['sequence', 'peptide', 'helm'];
const METHODS = ['kalign', 'muscle', 'clustal'];
const SEQUENCES = ['MKTAYIAKQR', 'MKTAYIAQR', 'MKAYIAKQRQ', 'MTAYIAKQR', 'MKTAYIKQRQ'];

function align(method: string, keepGaps: boolean): string[] {
  const width = Math.max(...SEQUENCES.map((s) => s.length));
  return SEQUENCES.map((s, i) => {
    const padded = s.padEnd(width, keepGaps ? '-' : ' ');
    return `${i + 1}  ${method === 'muscle' ? padded.split('').reverse().join('') : padded}`;
  });
}

export function msaWorkbenchPage(): HTMLElement {
  // --- toolbar ------------------------------------------------------------------------------------
  const settingsOn = signal(false);
  const toolbar = new Toolbar({ariaLabel: 'Workbench toolbar'});
  toolbar.name = 'workbenchToolbar';
  toolbar.addButton('Open', () => openDatasetDialog(), {icon: '📂', tooltip: 'Open a dataset'});
  toolbar.addButton('Run', () => openMsaDialog(), {icon: '▶', tooltip: 'Run the alignment'});
  const toolbarSave = button('Save', () => notify.info('Saved from the toolbar'));
  toolbarSave.dataset.u2Name = 'toolbarSave';
  toolbar.add(toolbarSave);
  toolbar.addSeparator();
  toolbar.addToggle('Settings', settingsOn, {tooltip: 'Show the settings panel'});

  // --- alignment form -----------------------------------------------------------------------------
  const nameInput = new TextInput({label: 'Name', name: 'alignmentName', placeholder: 'My alignment'});
  const sequenceColumn = new ChoiceInput({label: 'Sequence column', name: 'sequenceColumn', items: SEQUENCE_COLUMNS, value: 'sequence'});
  const method = new ChoiceInput({label: 'Method', name: 'method', items: METHODS, value: 'kalign'});
  const gapPenalty = new NumberInput({label: 'Gap penalty', name: 'gapPenalty', value: 1, min: 0, max: 10});
  const keepGaps = new BoolInput({label: 'Keep gaps', name: 'keepGaps', value: true});
  const form = new Form({layout: 'normal'});
  form.name = 'alignmentForm';
  for (const input of [nameInput, sequenceColumn, method, gapPenalty, keepGaps])
    form.add(input);
  const formSave = button('Save', () => notify.info('Saved from the form'));
  formSave.dataset.u2Name = 'formSave';
  const runMsa = button('Run MSA', () => openMsaDialog(), {primary: true});
  const alignmentPanel = divV([h2('Alignment'), form.root, divH([runMsa, formSave], 'u2demo-msa-actions')], 'u2demo-msa-card');
  alignmentPanel.dataset.u2Name = 'alignmentPanel';

  // --- results ------------------------------------------------------------------------------------
  const status = span('No alignment yet', 'u2demo-msa-status');
  status.dataset.u2Name = 'status';
  const alignedList = new VirtualList<string>({itemHeight: 24, render: (item) => span(item, 'u2demo-msa-sequence')});
  alignedList.name = 'alignedList';
  alignedList.root.classList.add('u2demo-msa-list');
  alignedList.setItems([]);
  const log = div([], 'u2demo-msa-log');
  log.dataset.u2Name = 'log';
  log.textContent = 'Nothing has run yet.';
  const resultTabs = new TabStrip({tabs: [
    {id: 'alignment', label: 'Alignment', content: alignedList.root},
    {id: 'log', label: 'Log', content: log},
  ]});
  resultTabs.name = 'resultTabs';
  const results = divV([h2('Results'), status, resultTabs.root], 'u2demo-msa-card');
  results.dataset.u2Name = 'results';

  // --- settings -----------------------------------------------------------------------------------
  const settings = new Section({title: 'Settings'});
  settings.name = 'settingsPanel';
  settings.body.append(new ChoiceInput({label: 'Theme', name: 'theme', items: ['light', 'dark'], value: 'light'}).root);
  toolbar.effect(() => settings.root.hidden = !settingsOn.value);

  // --- dialogs ------------------------------------------------------------------------------------
  function openMsaDialog(): void {
    const column = new ChoiceInput({label: 'Sequence column', name: 'dialogSequenceColumn', items: SEQUENCE_COLUMNS,
      value: sequenceColumn.value.peek() ?? 'sequence'});
    const dialogMethod = new ChoiceInput({label: 'Method', name: 'dialogMethod', items: METHODS,
      value: method.value.peek() ?? 'kalign'});
    const dialogForm = new Form({layout: 'normal'});
    dialogForm.add(column);
    dialogForm.add(dialogMethod);
    new Dialog('MSA', {name: 'msaDialog'})
      .add(dialogForm.root)
      .onOK(() => {
        const chosen = dialogMethod.value.peek() ?? 'kalign';
        const rows = align(chosen, keepGaps.value.peek());
        alignedList.setItems(rows);
        status.textContent = `Aligned ${rows.length} sequences with ${chosen}`;
        log.textContent = [`method: ${chosen}`, `column: ${column.value.peek()}`,
          `gap penalty: ${gapPenalty.value.peek()}`, `keep gaps: ${keepGaps.value.peek()}`,
          `name: ${nameInput.value.peek() || '(unnamed)'}`].join('\n');
        notify.info('Alignment finished');
      })
      .onCancel(() => undefined)
      .show();
  }

  function openDatasetDialog(): void {
    const dataset = new ChoiceInput({label: 'Dataset', name: 'dataset', items: ['spgi', 'demog', 'peptides'], value: 'spgi'});
    new Dialog('Open dataset', {name: 'openDialog'})
      .add(dataset.root)
      .onOK(() => notify.info(`Opened ${dataset.value.peek()}`))
      .onCancel(() => undefined)
      .show();
  }

  const header = divH([h1('MSA workbench'), toolbar.root], 'u2demo-msa-header');
  const layout = div([alignmentPanel, divV([results, settings.root])], 'u2demo-msa-layout');
  const root = divV([header, layout], 'u2demo-page u2demo-msa');
  root.dataset.u2Name = 'msaWorkbench';
  return root;
}
