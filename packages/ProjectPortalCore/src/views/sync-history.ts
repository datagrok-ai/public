/** A portfolio's sync history: its schedule, its owners, and each run with what it changed, newest
 * first. Opened from the board's Sync history button; `syncHistory(portfolio)` returns the runs as a
 * table. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {badge, button, divH, divV, span} from '@datagrok-libraries/u2';
import {PortfolioRow, SyncRunRow} from '../generated/db';
import {dueState, syncRunsOf} from '../sync-runs';
import {describeChange, describeSchedule, ProgramChange} from '../model/sync-runs';
import {autoAdvance, autoAdvanceOf, summarizeAutoAdvanceResult} from '../auto-advance';
import {describeAutoAdvanceConfig} from '../model/auto-advance';
import {stagesOf} from '../portal';
import {userNames} from '../users';

const named = <T extends HTMLElement>(el: T, name: string): T => {
  el.dataset.u2Name = name;
  return el;
};

const utc = (v: any): string => {
  const t = v === null || v === undefined ? NaN : typeof v === 'string' ? Date.parse(v) : Number(v.valueOf());
  return Number.isNaN(t) ? '' : `${new Date(t).toISOString().slice(0, 16).replace('T', ' ')} UTC`;
};

/** The runs of [portfolio] as a table, newest first. */
export async function syncHistoryFrame(portfolio: PortfolioRow, limit = 200): Promise<DG.DataFrame> {
  const runs = await syncRunsOf(portfolio, limit);
  const names = await userNames(runs.map((r) => r.run_by));
  const S = DG.COLUMN_TYPE.STRING;
  const I = DG.COLUMN_TYPE.INT;
  const time = (v: any) => v === null || v === undefined ? null : new Date(typeof v === 'string' ? v : Number(v.valueOf()));
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'started', runs.map((r) => time(r.started_on))),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'finished', runs.map((r) => time(r.finished_on))),
    DG.Column.fromList(S, 'trigger', runs.map((r) => r.trigger)),
    DG.Column.fromList(S, 'run by', runs.map((r) => r.run_by ? names.get(r.run_by) ?? r.run_by : null)),
    DG.Column.fromList(I, 'created', runs.map((r) => r.created ?? 0)),
    DG.Column.fromList(I, 'updated', runs.map((r) => r.updated ?? 0)),
    DG.Column.fromList(I, 'flagged', runs.map((r) => r.flagged ?? 0)),
    DG.Column.fromList(I, 'unresolved', runs.map((r) => r.unresolved ?? 0)),
    DG.Column.fromList(S, 'unresolved spellings', runs.map((r) => r.unresolved_spellings?.split('|').join(', ') ?? null)),
    DG.Column.fromList(S, 'error', runs.map((r) => r.error ?? null)),
    DG.Column.fromList(S, 'notified', runs.map((r) => r.notified ?? null)),
    DG.Column.fromList(S, 'changes', runs.map((r) => programChanges(r).map(describeChange).join('\n') || null)),
  ]);
  df.name = `${portfolio.code} · Sync history`;
  return df;
}

const programChanges = (r: SyncRunRow): ProgramChange[] => (r.changes as {programs?: ProgramChange[]} | undefined)?.programs ?? [];

/** The history as the app shows it. */
export async function syncHistoryPanel(portfolio: PortfolioRow, limit = 20): Promise<HTMLElement> {
  const state = await dueState(portfolio);
  const runs = await syncRunsOf(portfolio, limit);
  const names = await userNames(runs.map((r) => r.run_by));
  const schedule = state.schedule ?
    `Scheduled: ${describeSchedule(state.schedule)} · next ${state.due ? 'due now' : utc(state.next)}` :
    state.scheduleError ? `Schedule "${portfolio.schedule}" cannot be read: ${state.scheduleError}` :
      'No schedule: syncs only when asked';
  const owners = (portfolio.owners ?? '').split('|').filter((o) => o.length > 0);
  const advance = autoAdvanceOf(portfolio);
  const through = advance ? (await stagesOf(portfolio.stage_model_id)).find((s) => s.code === advance.through)?.name : undefined;
  const preview = named(divV([], 'pp-units'), 'auto-advance-preview');
  const previewButton = named(button('Preview advance by data', () => {
    preview.replaceChildren(span('Checking each program\'s next gate…', 'pp-muted'));
    autoAdvance(portfolio.code, true)
      .then((r) => preview.replaceChildren(...summarizeAutoAdvanceResult(r).split('\n').map((l) => span(l.trim(), 'pp-muted'))))
      .catch((e) => preview.replaceChildren(span(e?.message ?? String(e), 'pp-reason')));
  }), 'auto-advance-preview-button');
  const rows = runs.map((r, i) => {
    const changes = programChanges(r);
    const overlaps = (r.changes as {overlaps?: string[]} | undefined)?.overlaps ?? [];
    return named(divV([
      divH([
        span(utc(r.started_on)),
        badge(r.trigger === 'schedule' ? 'scheduled' : 'by hand', {variant: r.trigger === 'schedule' ? 'accent' : 'default'}),
        ...(r.run_by ? [span(names.get(r.run_by) ?? '', 'pp-muted')] : []),
        r.error ? badge('failed', {variant: 'error'}) : !r.finished_on ? badge('running or stopped', {variant: 'warning'}) :
          badge(`${r.created ?? 0} created · ${r.updated ?? 0} updated · ${r.flagged ?? 0} flagged`,
            {variant: (r.created ?? 0) + (r.updated ?? 0) + (r.flagged ?? 0) > 0 ? 'success' : 'default'}),
      ], 'pp-unit'),
      ...(r.error ? [span(r.error, 'pp-reason')] : []),
      ...(r.unresolved ? [span(`Unresolved stage spellings: ${(r.unresolved_spellings ?? '').split('|').join(', ')}`, 'pp-gap')] : []),
      ...(overlaps.length ? [span(`Overlaps: ${overlaps.join(', ')}`, 'pp-gap')] : []),
      ...changes.map((c) => span(describeChange(c), 'pp-muted')),
      ...(r.notified ? [span(`Owners: ${r.notified}`, 'pp-muted')] : []),
    ], 'pp-attention-item'), `sync-run-${i}`);
  });
  return named(divV([
    named(span(schedule), 'sync-schedule'),
    named(span(owners.length ? `Owners: ${owners.join(', ')}` : 'No owners declared', 'pp-muted'), 'sync-owners'),
    named(span(describeAutoAdvanceConfig(advance, through), 'pp-muted'), 'sync-auto-advance'),
    ...(advance ? [previewButton, preview] : []),
    span(`Runs (${runs.length}${runs.length === limit ? '+' : ''})`, 'pp-heading'),
    runs.length ? named(divV(rows, 'pp-units'), 'sync-runs') : named(span('Never synced', 'pp-muted'), 'no-sync-runs'),
  ], 'pp-sync-history'), 'sync-history');
}

/** Opens the history of [portfolio] in a dialog. */
export async function openSyncHistory(portfolio: PortfolioRow): Promise<void> {
  const panel = await syncHistoryPanel(portfolio);
  const dlg = ui.dialog({title: `${portfolio.name}: sync history`}).add(panel);
  dlg.root.dataset.u2Name = 'sync-history-dialog';
  dlg.show({width: 560});
}
