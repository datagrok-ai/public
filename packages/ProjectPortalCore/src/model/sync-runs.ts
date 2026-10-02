/** Scheduled sync (F-01): when a portfolio's sync is due, and what one run changed. Pure: the
 * schedule is read against the last run's start, and a run's plans are folded into the counts the
 * sync history shows and the per-program changes later steps (data-driven stage advance, new-data
 * notification) build on. Times are UTC. */
import type {SyncPlan} from './sync-plan';
import type {WorkUnitPlan} from './work-units';
import type {CompoundPlan} from './compounds';

const MINUTE = 60 * 1000;
const UNITS: Record<string, number> = {m: MINUTE, h: 60 * MINUTE, d: 24 * 60 * MINUTE};
/** The shortest interval a schedule may ask for: a sync reads every feed of the portfolio. */
export const MIN_INTERVAL_MS = 15 * MINUTE;

const ALIASES: Record<string, string> = {'@hourly': '0 * * * *', '@daily': '0 0 * * *', '@weekly': '0 0 * * 0'};
const CRON_FIELDS: [string, number, number][] = [['minute', 0, 59], ['hour', 0, 23], ['day of month', 1, 31],
  ['month', 1, 12], ['day of week', 0, 7]];

export type Schedule =
  {kind: 'interval', ms: number, text: string} |
  {kind: 'cron', fields: Set<number>[], restricted: [boolean, boolean], text: string};

function cronField(raw: string, [name, lo, hi]: [string, number, number]): Set<number> {
  const out = new Set<number>();
  for (const part of raw.split(',')) {
    const m = /^(\*|(\d+)(?:-(\d+))?)(?:\/(\d+))?$/.exec(part);
    if (!m)
      throw new Error(`${name} "${part}" is not a number, range, list or step`);
    const from = m[1] === '*' ? lo : Number(m[2]);
    const to = m[1] === '*' ? hi : m[3] !== undefined ? Number(m[3]) : m[4] !== undefined ? hi : from;
    const step = m[4] !== undefined ? Number(m[4]) : 1;
    if (from < lo || to > hi || from > to || step < 1)
      throw new Error(`${name} "${part}" is outside ${lo}-${hi}`);
    for (let v = from; v <= to; v += step)
      out.add(name === 'day of week' && v === 7 ? 0 : v);
  }
  return out;
}

/** A schedule as a portfolio file writes it: an interval (`every 6h`, `30m`, `1d`) or a five-field
 * cron line in UTC (`0 6 * * 1-5`, `@daily`). Throws with the reason when it cannot be read. */
export function parseSchedule(text: string): Schedule {
  const s = String(text ?? '').trim();
  const every = /^(?:every\s+)?(\d+)\s*(m|min|mins|minutes?|h|hours?|d|days?)$/i.exec(s);
  if (every) {
    const ms = Number(every[1]) * UNITS[every[2][0].toLowerCase()];
    if (ms < MIN_INTERVAL_MS)
      throw new Error(`"${s}" is shorter than ${MIN_INTERVAL_MS / MINUTE} minutes`);
    return {kind: 'interval', ms, text: s};
  }
  const parts = (ALIASES[s.toLowerCase()] ?? s).split(/\s+/);
  if (parts.length !== 5)
    throw new Error(`"${s}" is neither an interval (every 6h) nor a five-field cron line (0 6 * * 1-5)`);
  const fields = parts.map((p, i) => cronField(p, CRON_FIELDS[i]));
  return {kind: 'cron', fields, restricted: [parts[2] !== '*', parts[4] !== '*'], text: s};
}

/** Every problem with a portfolio's schedule and owners, as sentences. */
export function validateSchedule(schedule: string | undefined, owners: string[] | undefined): string[] {
  const errors: string[] = [];
  if (schedule !== undefined) {
    try {
      parseSchedule(schedule);
    } catch (e: any) {
      errors.push(`schedule: ${e.message}`);
    }
  }
  for (const o of owners ?? []) {
    if (!o || typeof o !== 'string' || o.includes('|'))
      errors.push(`owner "${o}" is not a user login`);
  }
  return errors;
}

/** The first minute strictly after [after] that a cron schedule names. */
function nextCron(s: Extract<Schedule, {kind: 'cron'}>, after: number): number {
  const [minutes, hours, days, months, weekdays] = s.fields;
  const [domSet, dowSet] = s.restricted;
  let t = Math.floor(after / MINUTE) * MINUTE + MINUTE;
  // at most five years of skips: month, then day, then hour, then minute
  for (let guard = 0; guard < 200000; guard++) {
    const d = new Date(t);
    const y = d.getUTCFullYear();
    const mo = d.getUTCMonth();
    if (!months.has(mo + 1)) {
      t = Date.UTC(y, mo + 1, 1);
      continue;
    }
    const dom = days.has(d.getUTCDate());
    const dow = weekdays.has(d.getUTCDay());
    // cron: when both day fields are restricted, either one matching is enough
    const dayOk = domSet && dowSet ? dom || dow : domSet ? dom : dowSet ? dow : true;
    if (!dayOk) {
      t = Date.UTC(y, mo, d.getUTCDate() + 1);
      continue;
    }
    if (!hours.has(d.getUTCHours())) {
      t = Date.UTC(y, mo, d.getUTCDate(), d.getUTCHours() + 1);
      continue;
    }
    if (!minutes.has(d.getUTCMinutes())) {
      t += MINUTE;
      continue;
    }
    return t;
  }
  throw new Error(`schedule "${s.text}" names no time in the next five years`);
}

/** When the next run falls: after the last run's start, or, never run, now. */
export function nextRun(schedule: Schedule, lastStarted: number | null, now: number): number {
  if (lastStarted === null)
    return now;
  return schedule.kind === 'interval' ? lastStarted + schedule.ms : nextCron(schedule, lastStarted);
}

/** Whether a portfolio's sync is due at [now]: never run, or its next run has come. Any run counts,
 * by schedule or by hand, and a failed one too: a failure waits for the next slot, not the next tick. */
export function isDue(schedule: Schedule, lastStarted: number | null, now: number): boolean {
  return nextRun(schedule, lastStarted, now) <= now;
}

/** In words, for the sync history: `every 6h` / `cron 0 6 * * 1-5 (UTC)`. */
export function describeSchedule(s: Schedule): string {
  return s.kind === 'interval' ? (/^every/i.test(s.text) ? s.text : `every ${s.text}`) : `cron ${s.text} (UTC)`;
}

// ---------------------------------------------------------------------------------------------
// What a run changed

/** One program's changes in a run: the hook data-driven stage advance and new-data notification read. */
export interface ProgramChange {
  program: string;
  created?: boolean;
  /** Program fields the registry changed. */
  fields?: string[];
  /** Stage ids, from and to. */
  stage?: {from: string | null, to: string | null};
  /** No longer in the feed (flagged, kept). */
  missing?: boolean;
  units?: {created: string[], updated: string[], missing: string[]};
  compounds?: {added: string[], roleChanges: {compound: string, from: string | null, to: string}[], missing: string[]};
}

export interface SyncRunCounts {
  created: number;
  updated: number;
  /** Programs, work units and compound links flagged in this run as no longer in their feed. */
  flagged: number;
  unresolved: number;
  /** `CODE="spelling"`, one per unresolved stage spelling. */
  unresolvedSpellings: string[];
  /** `CODE also in OTHER`: programs another portfolio's feed claims too. Standing, like spellings. */
  overlaps: string[];
}

export interface RunInput {
  programs: SyncPlan[];
  units: WorkUnitPlan[];
  compounds: CompoundPlan | null;
  /** Program ids to codes, for work units inserted under a program. */
  programCodes: Map<string, string>;
}

/** The counts and per-program changes of one run. */
export function summarizeRun(input: RunInput): {counts: SyncRunCounts, changes: ProgramChange[]} {
  const by = new Map<string, ProgramChange>();
  const of = (program: string) => {
    const c = by.get(program) ?? {program};
    by.set(program, c);
    return c;
  };
  const unitsOf = (program: string) => of(program).units ??= {created: [], updated: [], missing: []};
  const compoundsOf = (program: string) => of(program).compounds ??= {added: [], roleChanges: [], missing: []};
  const split = (key: string) => {
    const i = key.indexOf('/');
    return i < 0 ? [key, ''] : [key.slice(0, i), key.slice(i + 1)];
  };
  const counts: SyncRunCounts = {created: 0, updated: 0, flagged: 0, unresolved: 0, unresolvedSpellings: [], overlaps: []};

  for (const p of input.programs) {
    for (const v of p.inserts)
      of(v.code).created = true;
    for (const u of p.updates)
      of(u.code).fields = Object.keys(u.changes).sort();
    for (const s of p.stageChanges)
      of(s.code).stage = {from: s.fromStageId, to: s.toStageId};
    for (const m of p.missing)
      of(m.code).missing = true;
    counts.created += p.inserts.length;
    counts.updated += p.updates.length;
    counts.flagged += p.missing.length;
    counts.overlaps.push(...p.overlaps.map((o) => `${o.code} also in ${o.otherPortfolio}`));
    counts.unresolved += p.unresolvedStages.length;
    counts.unresolvedSpellings.push(...p.unresolvedStages.map((u) => `${u.code}="${u.spelling}"`));
  }
  for (const w of input.units) {
    for (const v of w.inserts)
      unitsOf(input.programCodes.get(v.program_id) ?? '?').created.push(`${w.kind} ${v.code}`);
    for (const u of [...w.updates, ...w.returned]) {
      const [program, code] = split(u.key);
      unitsOf(program).updated.push(`${w.kind} ${code}`);
    }
    for (const m of w.missing) {
      const [program, code] = split(m.key);
      unitsOf(program).missing.push(`${w.kind} ${code}`);
    }
    counts.created += w.inserts.length;
    counts.updated += w.updates.length + w.returned.length;
    counts.flagged += w.missing.length;
  }
  const c = input.compounds;
  if (c) {
    for (const l of c.linkInserts)
      compoundsOf(l.programCode).added.push(l.compoundCode);
    for (const r of c.roleChanges)
      compoundsOf(r.programCode).roleChanges.push({compound: r.compoundCode, from: r.from, to: r.to});
    for (const m of c.linkMissing) {
      const [program, code] = split(m.key);
      compoundsOf(program).missing.push(code);
    }
    counts.created += c.compoundInserts.length + c.linkInserts.length;
    counts.updated += c.compoundUpdates.length + c.roleChanges.length + c.linkReturned.length;
    counts.flagged += c.linkMissing.length;
  }
  const changes = [...by.values()].sort((a, b) => a.program.localeCompare(b.program));
  return {counts, changes};
}

/** One line per program change, for the sync history and the owners' notification. */
export function describeChange(c: ProgramChange): string {
  const parts: string[] = [];
  if (c.created)
    parts.push('new');
  if (c.fields?.length)
    parts.push(`changed ${c.fields.join(', ')}`);
  if (c.stage)
    parts.push('stage moved');
  if (c.missing)
    parts.push('no longer in the feed');
  const u = c.units;
  if (u) {
    if (u.created.length)
      parts.push(`new ${u.created.join(', ')}`);
    if (u.updated.length)
      parts.push(`updated ${u.updated.join(', ')}`);
    if (u.missing.length)
      parts.push(`missing ${u.missing.join(', ')}`);
  }
  const k = c.compounds;
  if (k) {
    if (k.added.length)
      parts.push(`compounds added ${k.added.join(', ')}`);
    for (const r of k.roleChanges)
      parts.push(`${r.compound} ${r.from ?? '—'} → ${r.to}`);
    if (k.missing.length)
      parts.push(`compounds missing ${k.missing.join(', ')}`);
  }
  return `${c.program}: ${parts.join('; ')}`;
}

/** The run in one line: `3 created, 1 updated, 0 flagged, 1 unresolved spelling (EXS-XYZ9="Lead Gen")`. */
export function describeCounts(c: SyncRunCounts, error?: string | null): string {
  const line = `${c.created} created, ${c.updated} updated, ${c.flagged} flagged, ${c.unresolved} unresolved ` +
    `spelling${c.unresolved === 1 ? '' : 's'}${c.unresolvedSpellings.length ? ` (${c.unresolvedSpellings.join(', ')})` : ''}` +
    (c.overlaps.length ? `, ${c.overlaps.length} overlap${c.overlaps.length === 1 ? '' : 's'} (${c.overlaps.join(', ')})` : '');
  return error ? `failed: ${error}` : line;
}

/** What stays wrong from run to run: unresolved spellings and overlaps. */
export const standingIssues = (c: Pick<SyncRunCounts, 'unresolvedSpellings' | 'overlaps'>): string[] =>
  [...c.unresolvedSpellings, ...c.overlaps].sort();

/** Whether the owners hear about a run: it failed, changed or flagged something, or its standing
 * issues (unresolved spellings, overlaps) differ from the run before. A quiet run that finds what
 * the last one found says nothing. */
export function worthNotifying(run: {counts: SyncRunCounts, error?: string | null}, previousIssues: string[] | null): boolean {
  if (run.error)
    return true;
  const c = run.counts;
  if (c.created + c.updated + c.flagged > 0)
    return true;
  return standingIssues(c).join('|') !== [...(previousIssues ?? [])].sort().join('|');
}
