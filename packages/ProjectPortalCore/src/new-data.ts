/** New-data notification against the platform (F-11): follows in `follow`, one row per user and
 * program or work unit; each sync (and `checkNewData`) reads the program dataset of every followed
 * program once, compares it with each follow's snapshot, tells followers in the platform's
 * notifications, and moves the snapshots on. The rules live in ./model/new-data. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {FollowRow, FollowWorkUnitKind, portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {portfolioByCode, programData} from './portal';
import {compare, DataSnapshot, describeNews, followLabel, isNews, rowsOf} from './model/new-data';
import {DataRow, ProgramDataConfig} from './model/conditions';
import {notifyUser} from './notify';
import {filterValue, ROW_LIMIT} from './rows';

export interface WorkUnitRef {
  kind: FollowWorkUnitKind;
  code: string;
}

const me = () => DG.User.current().id;
const watermarkOf = (p: PortfolioRow) => (p.program_data as ProgramDataConfig | undefined)?.watermark ?? null;
const isUnit = (f: FollowRow, unit?: WorkUnitRef | null) =>
  unit ? f.work_unit_kind === unit.kind && f.work_unit_code === unit.code : !f.work_unit_code;
/** Only a series has rows of its own in the program dataset (its `series` column); a follow of any
 * other work unit hears about the program. */
const seriesOf = (f: Pick<FollowRow, 'work_unit_kind' | 'work_unit_code'>) =>
  f.work_unit_kind === 'series' ? f.work_unit_code ?? null : null;

const isActive = (f: FollowRow) => f.active !== false;

/** Every follow row of the current user for a program, unfollowed ones too. */
const myRows = (programId: string): Promise<FollowRow[]> =>
  portalDb.follows.query({filter: `program_id = ${filterValue(programId)} and user_id = ${filterValue(me())}`, limit: ROW_LIMIT});

/** The current user's follows of a program: the program itself and any of its work units. */
export async function myFollows(programId: string): Promise<FollowRow[]> {
  return (await myRows(programId)).filter(isActive);
}

/** Every active follow of the current user, across programs: what "my feed" reads. */
export async function myFollowsAll(): Promise<FollowRow[]> {
  return (await portalDb.follows.query({filter: `user_id = ${filterValue(me())}`, limit: ROW_LIMIT})).filter(isActive);
}

async function programByCode(portfolioCode: string, programCode: string): Promise<{portfolio: PortfolioRow, program: ProgramRow}> {
  const portfolio = await portfolioByCode(portfolioCode);
  const program = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(programCode)}`, limit: ROW_LIMIT}))[0];
  if (!program)
    throw new Error(`${portfolioCode} has no program ${programCode}, or it is not shared with you`);
  return {portfolio, program};
}

/** Follows a program, or one of its work units, as the current user. The data as it stands now is
 * the baseline, so what is already there is never news; when it cannot be read, the next check
 * takes the baseline. Following twice is a no-op. */
export async function follow(portfolio: PortfolioRow, program: ProgramRow, unit?: WorkUnitRef | null): Promise<FollowRow> {
  const existing = (await myRows(program.id)).find((f) => isUnit(f, unit));
  if (existing && isActive(existing))
    return existing;
  let base: DataSnapshot | null = null;
  try {
    const rows = await programData(portfolio, program.code);
    if (rows)
      base = compare(null, rowsOf(rows, unit?.kind === 'series' ? unit.code : null), watermarkOf(portfolio)).now;
  } catch (e: any) {
    grok.log.warning(`Project Portal: no baseline for following ${program.code}: ${e?.message ?? e}`);
  }
  const now = new Date().toISOString();
  const snapshots = {checked: base ?? undefined, checked_on: base ? now as any : undefined,
    seen: base ?? undefined, seen_on: base ? now as any : undefined};
  // an unfollowed row is turned back on, with a fresh baseline
  if (existing) {
    // and without the news told before the unfollow, which would read as unseen
    await portalDb.follows.update(existing.id, {active: true, ...snapshots, checked: base as any, seen: base as any,
      checked_on: now as any, seen_on: now as any, last_news: null as any, news_on: null as any});
    return portalDb.follows.get(existing.id);
  }
  const [r] = await portalDb.follows.insert({user_id: me(), program_id: program.id,
    work_unit_kind: unit?.kind, work_unit_code: unit?.code, active: true, ...snapshots});
  return portalDb.follows.get((r as any).existingId ?? r.id);
}

/** Stops following; returns how many follows were turned off. The row stays, inactive. */
export async function unfollow(program: ProgramRow, unit?: WorkUnitRef | null): Promise<number> {
  const mine = (await myFollows(program.id)).filter((f) => isUnit(f, unit));
  for (const f of mine)
    await portalDb.follows.update(f.id, {active: false});
  return mine.length;
}

export async function followByCode(portfolioCode: string, programCode: string, series?: string | null,
  stop = false): Promise<string> {
  const {portfolio, program} = await programByCode(portfolioCode, programCode);
  const unit: WorkUnitRef | null = series ? {kind: 'series', code: series} : null;
  const label = followLabel(program.code, unit ?? undefined);
  if (stop)
    return `${label}: ${(await unfollow(program, unit)) ? 'no longer followed' : 'was not followed'}`;
  const f = await follow(portfolio, program, unit);
  return `${label}: followed${f.checked ? ` (${(f.checked as DataSnapshot).rows} rows now)` : ''}`;
}

export interface NewDataCheck {
  /** Follows looked at. */
  follows: number;
  /** Per follow with news: who, and the line. */
  news: {user: string, program: string, line: string}[];
  /** Per user told: how it went. */
  notified: string[];
  /** Programs whose data could not be read. */
  errors: string[];
  summary: string;
}

/** Checks every followed program of a portfolio for new data, tells the followers, and moves each
 * follow's snapshot on. Programs nobody follows are not read. */
export async function checkNewDataOf(portfolio: PortfolioRow, opts: {notify?: boolean} = {}): Promise<NewDataCheck> {
  const programs = await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)}`, limit: ROW_LIMIT});
  const ids = programs.map((p) => `"${p.id}"`).join(', ');
  const follows = ids ? (await portalDb.follows.query({filter: `program_id in (${ids})`, limit: ROW_LIMIT})).filter(isActive) : [];
  const out: NewDataCheck = {follows: follows.length, news: [], notified: [], errors: [], summary: ''};
  if (!follows.length) {
    out.summary = `${portfolio.code} new data: nobody follows its programs`;
    return out;
  }
  const watermark = watermarkOf(portfolio);
  const now = new Date().toISOString();
  let baselines = 0;
  for (const program of programs) {
    const mine = follows.filter((f) => f.program_id === program.id);
    if (!mine.length)
      continue;
    let rows: DataRow[] | null;
    try {
      rows = await programData(portfolio, program.code);
    } catch (e: any) {
      out.errors.push(`${program.code} (${e?.message ?? e})`);
      continue;
    }
    if (!rows) {
      out.errors.push(`${program.code} (${portfolio.code} declares no program data)`);
      continue;
    }
    for (const f of mine) {
      const n = compare(f.checked as DataSnapshot | undefined, rowsOf(rows, seriesOf(f)), watermark);
      const line = describeNews(followLabel(program.code, {kind: f.work_unit_kind, code: f.work_unit_code}), n);
      if (n.baseline)
        baselines++;
      // one follow that cannot be written does not cost the others their news
      try {
        await portalDb.follows.update(f.id, {checked: n.now, checked_on: now as any,
          ...(line ? {last_news: line, news_on: now as any} : {})});
      } catch (e: any) {
        out.errors.push(`${program.code} (follow not updated: ${e?.message ?? e})`);
        continue;
      }
      if (line)
        out.news.push({user: f.user_id, program: program.code, line});
    }
  }
  const users = [...new Set(out.news.map((x) => x.user))];
  if (opts.notify !== false) {
    for (const u of users) {
      const lines = out.news.filter((x) => x.user === u).map((x) => x.line);
      const r = await notifyUser(u, [`New data in ${portfolio.name} (${portfolio.code})`, ...lines].join('\n'));
      const who = await grok.dapi.users.find(u).then((x) => x.login).catch(() => u);
      out.notified.push(`${who}: ${r.detail}`);
    }
  }
  const lines = out.news.map((x) => x.line);
  out.summary = `${portfolio.code} new data: ` + (lines.length ? [...new Set(lines)].join('; ') : 'nothing new') +
    ` (${follows.length} follow(s)${baselines ? `, ${baselines} baseline(s) taken` : ''})` +
    (out.notified.length ? `; followers: ${out.notified.join(', ')}` : '') +
    (out.errors.length ? `; not read: ${out.errors.join(', ')}` : '');
  return out;
}

/** What landed in a program's data since the current user's last visit, per follow, and the visit
 * recorded. Null when the user does not follow the program or any of its work units. */
export async function whatsNew(portfolio: PortfolioRow, program: ProgramRow, rows: DataRow[] | null,
  follows?: FollowRow[]): Promise<{since: string | null, lines: string[]} | null> {
  const mine = follows ?? await myFollows(program.id);
  if (!mine.length || !rows)
    return null;
  const watermark = watermarkOf(portfolio);
  const since = mine.map((f) => f.seen_on ? String(f.seen_on.toISOString?.() ?? f.seen_on) : null)
    .filter((x): x is string => !!x).sort()[0] ?? null;
  const lines: string[] = [];
  const now = new Date().toISOString();
  for (const f of mine) {
    const n = compare(f.seen as DataSnapshot | undefined, rowsOf(rows, seriesOf(f)), watermark);
    const line = describeNews(followLabel(program.code, {kind: f.work_unit_kind, code: f.work_unit_code}), n);
    if (line)
      lines.push(line);
    await portalDb.follows.update(f.id, {seen: n.now, seen_on: now as any})
      .catch((e) => grok.log.warning(`Project Portal: visit not recorded: ${e?.message ?? e}`));
  }
  return {since, lines};
}

export {isNews};
