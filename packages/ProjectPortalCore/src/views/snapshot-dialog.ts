/** One program snapshot, read-only (F-32): when and by whom it was taken, the stage and status, the
 * state of each criterion, the runs selected, and the table as frozen, which opens as a table view
 * with editing off. Whether the table file still holds what was frozen is checked on open. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {divV, span} from '@datagrok-libraries/u2';
import {portalDb, ProgramSnapshotRow} from '../generated/db';
import {frozenCriteria, frozenRuns, snapshotTable, summaryOf} from '../snapshots';
import {describeSnapshot, snapshotDate} from '../model/snapshots';
import {installStyles} from './styles';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const CHECK_TEXT = {
  'intact': 'the file holds what was frozen',
  'changed': 'the file was changed after the snapshot: its SHA-256 differs',
  'missing': 'the file is gone',
  'not kept': 'not kept',
} as const;

function table(head: string[], rows: string[][], name: string): HTMLElement {
  const t = named(document.createElement('table'), name);
  t.className = 'pp-card-table';
  const h = t.insertRow();
  for (const x of head)
    h.append(Object.assign(document.createElement('th'), {textContent: x}));
  for (const r of rows) {
    const tr = t.insertRow();
    for (const x of r)
      tr.insertCell().textContent = x;
  }
  return t;
}

/** Opens the frozen table as a table view the user cannot edit. */
export async function openSnapshotTable(s: ProgramSnapshotRow): Promise<DG.TableView | null> {
  const {df, check} = await snapshotTable(s);
  if (!df) {
    grok.shell.warning(`The snapshot's table: ${CHECK_TEXT[check]}`);
    return null;
  }
  const tv = grok.shell.addTableView(df);
  tv.grid.setOptions({allowEdit: false});
  if (check === 'changed')
    grok.shell.warning(`${df.name}: ${CHECK_TEXT.changed}`);
  return tv;
}

export async function openSnapshotDialog(snapshotId: string): Promise<DG.Dialog> {
  installStyles();
  const s = await portalDb.programSnapshots.get(snapshotId);
  const program = await portalDb.programs.get(s.program_id).catch(() => null);
  const by = s.taken_by ? (await grok.dapi.users.find(s.taken_by).catch(() => null))?.friendlyName ?? '' : '';
  const summary = summaryOf(s);
  const {check} = s.table_file ? await snapshotTable(s).catch(() => ({check: 'missing' as const})) : {check: 'not kept' as const};
  const criteria = frozenCriteria(s);
  const runs = frozenRuns(s);
  const open = named(ui.button('Open table', () => void openSnapshotTable(s).catch((e) => grok.shell.error(e?.message ?? String(e)))),
    'snapshot-open-table');
  open.title = 'The program dataset as frozen, in a table view you cannot edit';
  (open as HTMLButtonElement).disabled = !s.table_file || check === 'missing';
  const dlg = ui.dialog({title: `Snapshot of ${program?.code ?? 'a program'} · ${snapshotDate(summary.takenOn)}`});
  dlg.root.dataset.u2Name = 'snapshot-dialog';
  dlg.add(named(divV([
    named(span(describeSnapshot(summary)), 'snapshot-summary'),
    span(`${s.label}${by ? ` · taken by ${by}` : ''}`, 'pp-muted'),
    named(span(`Stage: ${s.stage_name ?? 'unstaged'}${s.stage_code ? ` (${s.stage_code})` : ''} · status ${String(s.status ?? '').replace(/_/g, ' ')}`),
      'snapshot-stage'),
    ...(s.note ? [named(span(s.note, 'pp-muted'), 'snapshot-note')] : []),
    span('Criteria', 'pp-heading'),
    criteria.length ? table(['Stage', 'Criterion', 'State', 'What the data said'],
      criteria.map((c) => [c.stage, `${c.name}${c.required ? '' : ' (optional)'}`, c.state.replace(/_/g, ' '),
        c.summary ?? (c.basis === 'evaluated' ? 'from data' : '')]), 'snapshot-criteria') :
      span('No criteria apply', 'pp-muted'),
    span('Selected runs', 'pp-heading'),
    runs.length ? table(['Criterion', 'Run', 'Model', 'Why'], runs.map((r) => [r.criterion, r.title, r.model, r.rationale ?? '']),
      'snapshot-runs') : span('None', 'pp-muted'),
    span(s.dataset ? `Table: ${s.dataset}` : 'Table', 'pp-heading'),
    named(span(s.table_rows === undefined || s.table_rows === null ? 'No table' :
      `${s.table_rows} rows${s.table_file ? ` · ${s.table_file}` : ''} · ${CHECK_TEXT[check]}`, check === 'changed' ? 'pp-reason' : ''),
    'snapshot-table'),
    open,
  ], 'pp-units'), 'snapshot-content'));
  dlg.show({width: 600});
  return dlg;
}
