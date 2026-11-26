/* The dialog that plans a tumor model the atlas has never seen from its own history: the lab's or the provider's
   pilot, or earlier studies, read from a file (the study reader: long or wide tables, Excel, the volume formula) or
   from a table open in Datagrok; its control group, a readable name, four short questions and, optionally, the
   pilot's take. Registers the model (history-models.ts) and gives its key. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {controlByName} from './arm-name';
import {SavedHistory, historyFit, registerHistory} from './history-models';
import type {Answer, HistoryRow} from './model-history';
import {fileText, readFile} from './open-study';
import {previousOfTable} from './previous-studies';
import {dialogLine, okNamed, showCentred, modal} from './ui';

/** A measurement as the readers give it. */
type ReadRow = {animal: string, group: string | null, day: number | null, volume: number | null};

/** The control group's measurements with a day and a volume above 0; with no group column, every row. */
export const controlRows = (rows: ReadRow[], control: string | null): HistoryRow[] => rows
  .filter((r) => (control === null || r.group === control) && r.day !== null && r.volume !== null && r.volume > 0)
  .map((r) => ({animal: r.animal, day: r.day!, volume: r.volume!}));

/** A model's name from its file or table: without the extension or a trailing «history». */
export const defaultName = (from: string) =>
  from.replace(/\.(csv|tsv|txt|xlsx)$/i, '').replace(/[\s_-]*history$/i, '').trim() || from;

/** Whether the pilot's take, percent, is one (empty: none given). */
export const takeOk = (take: number | null) => take === null || (take >= 1 && take <= 100);

const ERROR_COLOR = 'var(--red-3, #c62828)';

/** The groups of the rows, in their order. */
const groupsOf = (rows: ReadRow[]) => [...new Set(rows.map((r) => r.group).filter((g): g is string => g !== null))];

const ANSWERS: {[text: string]: Answer} = {'Yes': 'yes', 'No': 'no', 'Not sure': 'unknown'};
const PASSAGES: {[text: string]: SavedHistory['source']['passage']} =
  {'P0–P1': 'early', 'P2 or later': 'later', 'Not sure': 'unknown'};

/** Opens the dialog; `done` gets the registered model's key. `first`: a history already read (a table opened in
 *  Plan's Open a study), shown at once. */
export function historyDialog(done: (key: string) => void, first?: {rows: ReadRow[], from: string}): void {
  let read: {rows: ReadRow[], from: string} | null = null;
  // the name the dialog gave from the last file or table read ('' none yet)
  let named = '';
  const status = dialogLine('No history read yet.');
  const name = ui.input.string('Model name', {value: '', placeholder: 'e.g. PDX-LU-117',
    tooltipText: 'How lists, plans and reports name this model.'});
  const control = ui.input.choice<string>('Control group', {items: [],
    tooltipText: 'Only this group\'s growth is used.'});
  const lab = ui.input.choice('Run in your lab', {value: 'Not sure', items: Object.keys(ANSWERS),
    tooltipText: 'A history from another lab, or not known, gives wider ranges.'});
  const host = ui.input.choice('Same mouse strain and implant site', {value: 'Not sure',
    items: Object.keys(ANSWERS), tooltipText: 'Another strain or site gives wider ranges.'});
  const passage = ui.input.choice('Passage', {value: 'Not sure', items: Object.keys(PASSAGES),
    tooltipText: 'The tumor usually grows faster after its first passages.'});
  const frozen = ui.input.bool('From frozen material', {value: false,
    tooltipText: 'Tumors from frozen material may grow slower at first.'});
  const means = ui.input.bool('Rows are group means', {value: false,
    tooltipText: 'The history gives a group\'s mean per day, not each animal.'});
  const perGroup = ui.input.int('Animals per group', {nullable: true, min: 1,
    tooltipText: 'For group means: how many animals each mean is of. Empty: PMX Studio assumes 3.'});
  const take = ui.input.float('Tumors that grew, %', {nullable: true, min: 1, max: 100,
    tooltipText: 'Optional, from the pilot: of the tumors implanted, the share that grew. It sets only how many ' +
      'animals to implant. Empty: PMX Studio\'s PDX reference.'});
  // (the input's min and max only bound its spinner: a typed 250 stays, so the field and OK check it)
  take.addValidator((v) => takeOk(v?.trim() ? Number(v) : null) ? null :
    'A share of the tumors implanted: from 1 to 100%');
  const tables = grok.shell.tables.filter((t) => previousOfTable(t) !== null);
  const table = ui.input.choice<string>('Or a table open in Datagrok', {nullable: true,
    items: tables.map((t) => t.name)});

  const d = modal('history', 'Model not in the atlas');
  if (!d)
    return;
  let ok: HTMLButtonElement | null = null; // named once onOK has made it
  // the control's measurements of the history read (with no group column, every row); none read: none
  const controls = () => read ? controlRows(read.rows, groupsOf(read.rows).length ? control.value : null) : [];
  const refresh = () => {
    const rows = controls();
    const animals = new Set(rows.map((r) => r.animal)).size;
    const days = rows.map((r) => r.day);
    if (read) {
      status.textContent = rows.length ? `${read.from}: ${animals} control ${animals === 1 ? 'curve' : 'curves'}, ` +
        `days ${Math.min(...days)}–${Math.max(...days)}.` : `${read.from}: no control measurements read.`;
    }
    perGroup.root.style.display = means.value ? '' : 'none';
    if (ok)
      ok.disabled = !rows.length || !name.value.trim() || !takeOk(take.value);
  };
  const show = (rows: ReadRow[], from: string) => {
    read = {rows, from};
    status.style.color = '';
    const groups = groupsOf(rows);
    control.items = groups;
    control.value = controlByName(groups) ?? groups[0] ?? null;
    control.root.style.display = groups.length > 1 ? '' : 'none';
    // the name follows the file or table read until the user types one of their own
    if (!name.value || name.value === named)
      name.value = named = defaultName(from);
    // one curve in the control: a group's mean, not an animal
    means.value = new Set(controls().map((r) => r.animal)).size === 1;
    refresh();
  };

  const file = document.createElement('input');
  file.type = 'file';
  file.accept = '.csv,.tsv,.txt,.xlsx';
  file.addEventListener('change', async () => {
    const f = file.files?.[0];
    if (!f)
      return;
    try {
      const t = await fileText(f);
      const r = readFile(f.name, t.text, {}, t.merged);
      if (r.kind !== 'table')
        throw new Error(r.kind === 'bad' ? r.why : 'Open a table of measurements, not a plan or a study.');
      if (r.read.error || r.read.questions.length) {
        throw new Error(r.read.error ?? 'PMX Studio needs to ask about this table: open it first with Open a ' +
          'study, then pick it here among the open tables.');
      }
      show(r.read.rows, f.name);
    } catch (e: any) {
      // nothing of the history read before stays: its groups, the name it gave and Use this model
      read = null;
      control.items = [];
      control.root.style.display = 'none';
      if (name.value === named)
        name.value = named = '';
      refresh();
      status.textContent = e?.message ?? 'The file could not be read.';
      status.style.color = ERROR_COLOR;
    }
  });
  const choose = ui.button('Choose a file…', () => file.click(), 'A table of tumor volumes over days: CSV, TSV ' +
    'or Excel, long or wide.');
  table.onChanged.subscribe(() => {
    const df = tables.find((t) => t.name === table.value);
    const p = df ? previousOfTable(df) : null;
    if (p)
      show(p.rows, df!.name);
  });
  name.onChanged.subscribe(() => refresh());
  control.onChanged.subscribe(() => refresh());
  means.onChanged.subscribe(() => refresh());
  take.onChanged.subscribe(() => refresh());

  d.add(dialogLine('Plan with this model\'s own growth: its implantation pilot, the provider\'s, or earlier ' +
    'studies.', 'PMX Studio weighs the history against its PDX reference: the more and closer the history, the ' +
    'more it counts.'))
    .add(ui.divH([choose, status], {style: {alignItems: 'center', gap: '8px'}}))
    .add(ui.inputs([...tables.length ? [table] : [], name, control, lab, host, passage, frozen, means, perGroup,
      take]));
  d.onOK(async () => {
    const h: SavedHistory = {name: name.value.trim(), take: take.value ?? null, rows: controls(),
      source: {sameLab: ANSWERS[lab.value!], sameHost: ANSWERS[host.value!], passage: PASSAGES[passage.value!],
        frozen: frozen.value, meansOnly: means.value, animalsPerGroup: means.value ? perGroup.value ?? null : null}};
    const key = registerHistory(h);
    try {
      await historyFit(key);
      done(key);
    } catch (e: any) {
      grok.shell.error(e?.message ?? 'The history could not be used.');
    }
  });
  ok = okNamed(d, 'Use this model');
  control.root.style.display = 'none';
  refresh();
  if (first)
    show(first.rows, first.from);
  showCentred(d, {width: 460});
}
