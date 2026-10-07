/** Scheduled sync against the platform (F-01): every sync of a portfolio, by schedule or by hand,
 * is recorded in `sync_log` with what it changed (linked into the portfolio's root, so only its
 * groups read it); `syncDuePortfolios` syncs only the portfolios
 * whose schedule is due; owners hear about a run that changed something in the platform's own
 * notifications. The rules live in ./model/sync-runs. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, SyncLogRow} from './generated/db';
import {portfolioByCode, summarizePlan, syncCompounds, syncPortfolio, syncWorkUnits} from './portal';
import {summarizeWorkUnitPlan} from './model/work-units';
import {summarizeCompoundPlan} from './model/compounds';
import {
  describeChange, describeCounts, describeSchedule, isDue, nextRun, parseSchedule, ProgramChange, Schedule,
  standingIssues, summarizeRun, SyncRunCounts, withAdvances, worthNotifying,
} from './model/sync-runs';
import {autoAdvance, autoAdvanceOf} from './auto-advance';
import {notifyUser} from './notify';
import {checkNewDataOf} from './new-data';
import {checkEvidenceChangesOf} from './evidence-tables';
import {myPortfolioIds, placeSyncLog} from './secured-rows';
import {filterValue, ROW_LIMIT} from './rows';

export type SyncTrigger = 'schedule' | 'manual';

export interface SyncRunResult {
  portfolio: string;
  trigger: SyncTrigger;
  /** The `sync_log` row; null when it could not be written (see `historyError`). */
  runId: string | null;
  historyError?: string;
  /** As syncPortfolio has always returned it: one line per feed. */
  report: string;
  counts: SyncRunCounts;
  /** Per program: what this run changed. The hook for data-driven stage advance and new-data notification. */
  changes: ProgramChange[];
  error: string | null;
  notified?: string;
  /** What followers were told after the run, one line per program or series, for the activity feed. */
  newData?: {program: string, line: string}[];
}

const ms = (v: any): number | null => {
  if (v === null || v === undefined || v === '')
    return null;
  const n = typeof v === 'string' ? Date.parse(v) : Number(v.valueOf());
  return Number.isNaN(n) ? null : n;
};

/** A portfolio's runs, newest first. */
export async function syncRunsOf(portfolio: PortfolioRow, limit = 50): Promise<SyncLogRow[]> {
  const mine = await myPortfolioIds();
  if (mine && !mine.has(portfolio.id))
    return [];
  const rows = await portalDb.syncLogs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, sort: '!started_on', limit});
  // the server's order is asked for, and kept whatever it returns
  return rows.sort((a, b) => (ms(b.started_on) ?? 0) - (ms(a.started_on) ?? 0));
}

const issuesOf = (run: SyncLogRow | undefined): string[] | null => !run ? null :
  [...(run.unresolved_spellings ?? '').split('|').filter((s) => s.length > 0),
    ...((run.changes as {overlaps?: string[]} | undefined)?.overlaps ?? [])].sort();

/** Syncs one portfolio (programs, then work units, then compounds) and records the run. The sync
 * never fails because its history could not be written; that is reported on the result. */
export async function recordedSync(code: string, trigger: SyncTrigger): Promise<SyncRunResult> {
  const portfolio = await portfolioByCode(code);
  const started = new Date().toISOString();
  let runId: string | null = null;
  let historyError: string | undefined;
  try {
    const [r] = await portalDb.syncLogs.insert({portfolio_id: portfolio.id, started_on: started, trigger,
      run_by: DG.User.current().id});
    runId = (r as any).existingId ?? r.id;
  } catch (e: any) {
    historyError = `sync history not recorded: ${e?.message ?? e}`;
  }
  // the row reaches the portfolio's groups through its root, and nobody else
  const unplaced = runId ? await placeSyncLog(runId, portfolio) : null;
  if (unplaced)
    historyError = `sync history not shared: ${unplaced}`;

  const empty: SyncRunCounts = {created: 0, updated: 0, flagged: 0, unresolved: 0, unresolvedSpellings: [], overlaps: []};
  let result: SyncRunResult = {portfolio: code, trigger, runId, historyError, report: '', counts: empty, changes: [], error: null};
  try {
    const programs = await syncPortfolio(code);
    const units = await syncWorkUnits(code);
    const compounds = await syncCompounds(code);
    const rows = await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT});
    const {counts, changes} = summarizeRun({programs, units, compounds, programCodes: new Map(rows.map((p) => [p.id, p.code]))});
    const report = [...programs.map(summarizePlan), ...units.map(summarizeWorkUnitPlan),
      ...(compounds ? [summarizeCompoundPlan(compounds)] : [])].join('\n');
    result = {...result, report, counts, changes};
  } catch (e: any) {
    result.error = e?.message ?? String(e);
  }
  // a portfolio that opts in lets its programs advance by data once the feeds are in
  if (!result.error && autoAdvanceOf(portfolio)) {
    try {
      const a = await autoAdvance(code, false);
      const moved = a.outcomes.filter((o) => o.to && !o.proposal);
      const proposed = a.outcomes.filter((o) => o.to && o.proposal);
      result.changes = withAdvances(result.changes, moved.map((o) => ({program: o.program, from: o.from?.code ?? null, to: o.to!.code})));
      const failed = a.outcomes.filter((o) => o.error);
      result.report += `\n${code} auto-advance: ${moved.length ? moved.map((o) => `${o.program} ${o.from?.code ?? '—'} → ${o.to!.code}`).join(', ') :
        proposed.length ? 'none moved' : `none of ${a.outcomes.length} program(s) met its next gate by data`}` +
        (failed.length ? `; not recorded: ${failed.map((o) => `${o.program} (${o.error})`).join(', ')}` : '') +
        (proposed.length ? `; waiting for approval: ${proposed.map((o) => `${o.program} ${o.from?.code ?? '—'} → ${o.to!.code}` +
          `${o.proposal === 'pending' ? ' (proposed earlier)' : ''}`).join(', ')}` : '');
    } catch (e: any) {
      result.report += `\n${code} auto-advance not applied: ${e?.message ?? e}`;
    }
  }

  // followers of a program or series hear about results that landed since the last check
  if (!result.error) {
    try {
      const n = await checkNewDataOf(portfolio);
      if (n.follows > 0)
        result.report += `\n${n.summary}`;
      const seen = new Set<string>();
      result.newData = n.news.filter((x) => !seen.has(x.line) && !!seen.add(x.line)).map((x) => ({program: x.program, line: x.line}));
    } catch (e: any) {
      grok.log.warning(`Project Portal: new data of ${code} not checked: ${e?.message ?? e}`);
    }
  }

  // generated documents whose evidence changed: marked stale, their author and followers told (F-41)
  if (!result.error) {
    try {
      const ev = await checkEvidenceChangesOf(portfolio);
      if (ev.documents > 0)
        result.report += `\n${ev.summary}`;
    } catch (e: any) {
      grok.log.warning(`Project Portal: evidence of ${code} not checked: ${e?.message ?? e}`);
    }
  }

  if (runId) {
    try {
      await portalDb.syncLogs.update(runId, {finished_on: new Date().toISOString() as any,
        created: result.counts.created, updated: result.counts.updated, flagged: result.counts.flagged,
        unresolved: result.counts.unresolved, unresolved_spellings: result.counts.unresolvedSpellings.join('|') || undefined,
        error: result.error ?? undefined, summary: result.report || undefined,
        changes: {programs: result.changes, overlaps: result.counts.overlaps,
          ...(result.newData?.length ? {newData: result.newData} : {})}});
    } catch (e: any) {
      result.historyError = `sync history not completed: ${e?.message ?? e}`;
    }
  }
  return result;
}

/** One run as owners read it: the counts, then a line per program it changed. */
export function runMessage(portfolio: {code: string, name: string}, r: Pick<SyncRunResult, 'counts' | 'changes' | 'error'>): string {
  const lines = [`${portfolio.name} (${portfolio.code}) sync: ${describeCounts(r.counts, r.error)}`];
  lines.push(...r.changes.slice(0, 20).map(describeChange));
  if (r.changes.length > 20)
    lines.push(`… and ${r.changes.length - 20} more programs`);
  return lines.join('\n');
}

/** Tells the portfolio's owners about a run, in the platform's notifications (the bell). The server
 * lets a user notify only themselves unless they may edit users, so a run by anyone else reaches the
 * owners it can and reports the rest. Never throws. */
export async function notifyOwners(portfolio: PortfolioRow, text: string): Promise<string> {
  const logins = (portfolio.owners ?? '').split('|').filter((o) => o.length > 0);
  if (logins.length === 0)
    return 'no owners declared';
  const out: string[] = [];
  for (const login of logins) {
    try {
      const user = await grok.dapi.users.filter(`login = "${login}"`).first().catch(() => null);
      if (!user || user.login !== login) {
        out.push(`${login}: no such user`);
        continue;
      }
      out.push(`${login}: ${(await notifyUser(user.id, text)).detail}`);
    } catch (e: any) {
      out.push(`${login}: not notified (${e?.message ?? e})`);
    }
  }
  return out.join('; ');
}

export interface DueState {
  portfolio: PortfolioRow;
  schedule: Schedule | null;
  /** Why a schedule could not be read. */
  scheduleError?: string;
  last: SyncLogRow | undefined;
  next: number | null;
  due: boolean;
}

/** Where a portfolio's schedule stands at [now]. */
export async function dueState(portfolio: PortfolioRow, now = Date.now()): Promise<DueState> {
  const last = (await syncRunsOf(portfolio, 1))[0];
  if (!portfolio.schedule)
    return {portfolio, schedule: null, last, next: null, due: false};
  try {
    const schedule = parseSchedule(portfolio.schedule);
    const started = ms(last?.started_on);
    return {portfolio, schedule, last, next: nextRun(schedule, started, now), due: isDue(schedule, started, now)};
  } catch (e: any) {
    return {portfolio, schedule: null, scheduleError: e?.message ?? String(e), last, next: null, due: false};
  }
}

/** Syncs the portfolios whose schedule is due at [now], discovery before development so lineage
 * resolves in one pass. Only portfolios that declare a schedule are looked at; [codes] narrows them
 * further. Owners hear about a run worth telling (see worthNotifying). */
export async function syncDue(codes?: string[], now = Date.now()): Promise<{ran: SyncRunResult[], skipped: string[]}> {
  const all = await portalDb.portfolios.query({limit: ROW_LIMIT});
  const wanted = codes?.length ? new Set(codes) : null;
  const order = (s?: string) => s === 'discovery' ? 0 : s === 'development' ? 1 : 2;
  const scheduled = all.filter((p) => p.schedule && (!wanted || wanted.has(p.code)))
    .sort((a, b) => order(a.setting) - order(b.setting) || a.code.localeCompare(b.code));
  const ran: SyncRunResult[] = [];
  const skipped: string[] = [];
  for (const code of wanted ?? []) {
    if (!all.some((p) => p.code === code))
      skipped.push(`${code}: not loaded`);
    else if (!scheduled.some((p) => p.code === code))
      skipped.push(`${code}: no schedule`);
  }
  for (const p of scheduled) {
    const state = await dueState(p, now);
    if (state.scheduleError) {
      skipped.push(`${p.code}: schedule "${p.schedule}" cannot be read: ${state.scheduleError}`);
      continue;
    }
    if (!state.due) {
      skipped.push(`${p.code}: next due ${new Date(state.next!).toISOString().slice(0, 16).replace('T', ' ')} UTC`);
      continue;
    }
    const r = await recordedSync(p.code, 'schedule');
    if (worthNotifying(r, issuesOf(state.last))) {
      r.notified = await notifyOwners(p, runMessage(p, r));
      if (r.runId)
        await portalDb.syncLogs.update(r.runId, {notified: r.notified}).catch(() => undefined);
    } else
      r.notified = 'nothing new to tell';
    ran.push(r);
  }
  return {ran, skipped};
}

/** The text syncDuePortfolios returns: a line per run, then what was skipped and why. */
export function summarizeDue(r: {ran: SyncRunResult[], skipped: string[]}): string {
  const lines = r.ran.map((x) => `${x.portfolio}: ${describeCounts(x.counts, x.error)}` +
    (x.notified ? `; owners: ${x.notified}` : '') + (x.historyError ? `; ${x.historyError}` : ''));
  if (r.skipped.length)
    lines.push(`not run: ${r.skipped.join('; ')}`);
  return lines.join('\n') || 'no portfolio declares a schedule';
}

export {describeSchedule, standingIssues};
