/** DMTA status (F-24, FR-105), discovery: where each compound of a program stands in the
 * design → make → test → analyze loop, designed → requested → made → tested → analyzed, and the
 * cycle time per series. Request objects (F-49) do not exist yet, so the status is derived from
 * what the portal has, one source per stage, and each source is pluggable:
 * - designed, requested: design sets (F-23), a proposal and its choice for synthesis;
 * - made: the compound registry projection (the compound is registered in the program);
 * - tested: a program dataset with a compound and a date column (the assay results);
 * - analyzed: an observation on the compound (F-26), or a role set by hand;
 * - anything else: functions the portfolio names (`dmta.sources`), e.g. a request system's status.
 * Pure: the caller reads the sources into events. */

export type DmtaStage = 'designed' | 'requested' | 'made' | 'tested' | 'analyzed';
export const DMTA_STAGES: DmtaStage[] = ['designed', 'requested', 'made', 'tested', 'analyzed'];

/** One fact: [compound] reached [stage], on a date when the source knows it. */
export interface DmtaEvent {
  /** The registered compound's code; a proposal not registered yet goes by its design label. */
  compound: string;
  stage: DmtaStage;
  /** ISO time; null when the source says it happened but not when (a registered compound). */
  on: string | null;
  /** Where it comes from: `design set KIN-A-2026-04-20`, `registry`, `Assay results`. */
  source: string;
  series?: string | null;
  /** The design label, when the compound was proposed in a design set. */
  label?: string | null;
}

/** A portfolio's DMTA settings (optional). */
export interface DmtaConfig {
  /** The dataset whose rows say a compound was tested, and on which date. Default: the first
   * program dataset with a compound column and a date column (run_date, tested_on, date). */
  tested?: {dataset: string, compoundColumn?: string, dateColumn?: string};
  /** Functions (Package:Function) that return more events: a table with compound, stage, on and,
   * optionally, source and series. Inputs `portfolio` and `program` (codes) when declared. */
  sources?: string[];
}

export function validateDmta(c: DmtaConfig | undefined, datasets: string[]): string[] {
  if (c === undefined)
    return [];
  const errors: string[] = [];
  if (c.tested !== undefined) {
    if (!c.tested?.dataset?.trim())
      errors.push('dmta: tested needs a dataset code');
    else if (!datasets.includes(c.tested.dataset))
      errors.push(`dmta: tested names dataset "${c.tested.dataset}", which the portfolio does not declare`);
  }
  for (const f of c.sources ?? []) {
    if (typeof f !== 'string' || !/^[^:\s]+:[^:\s]+$/.test(f))
      errors.push(`dmta: source "${f}" must be a function as Package:Function`);
  }
  return errors;
}

export const COMPOUND_COLUMNS = ['compound', 'compound_code', 'compound_id', 'code'];
export const DATE_COLUMNS = ['run_date', 'tested_on', 'test_date', 'date', 'measured_on'];

/** The compound and date columns of a dataset, by the names given or the usual names; null when
 * either is missing. */
export function testedColumns(columns: string[], asked?: {compoundColumn?: string, dateColumn?: string}):
  {compound: string, date: string} | null {
  const find = (names: string[]) => columns.find((c) => names.includes(c.toLowerCase())) ?? null;
  const compound = asked?.compoundColumn ? columns.find((c) => c === asked.compoundColumn) ?? null : find(COMPOUND_COLUMNS);
  const date = asked?.dateColumn ? columns.find((c) => c === asked.dateColumn) ?? null : find(DATE_COLUMNS);
  return compound && date ? {compound, date} : null;
}

const isoOf = (v: unknown): string | null => {
  if (v === null || v === undefined || v === '')
    return null;
  const d = typeof (v as any)?.toDate === 'function' ? (v as any).toDate() : v instanceof Date ? v : new Date(String(v));
  return Number.isNaN(d.getTime()) ? null : d.toISOString();
};

/** One `tested` event per compound: its first dated row. */
export function testedEvents(rows: Record<string, unknown>[], cols: {compound: string, date: string}, source: string): DmtaEvent[] {
  const first = new Map<string, string | null>();
  for (const r of rows) {
    const c = r[cols.compound];
    if (c === null || c === undefined || String(c).trim() === '')
      continue;
    const code = String(c).trim();
    const on = isoOf(r[cols.date]);
    const was = first.get(code);
    if (!first.has(code) || (on && (!was || on < was)))
      first.set(code, on);
  }
  return [...first.entries()].map(([compound, on]) => ({compound, stage: 'tested' as const, on, source}));
}

/** The compound a proposal became: its recorded code, else a registered compound with the same
 * SMILES (as written; no canonicalization). */
export function proposalCompound(p: {compound_code?: string | null, structure: string},
  compounds: {code: string, smiles?: string | null}[]): string | null {
  if (p.compound_code)
    return p.compound_code;
  const s = p.structure.trim();
  return compounds.find((c) => c.smiles && c.smiles.trim() === s)?.code ?? null;
}

/** One compound's loop: the date (or reached without a date) of each stage, and the furthest. */
export interface DmtaRow {
  compound: string;
  label: string | null;
  series: string | null;
  /** Per stage: the earliest date, '' when reached without a date, absent when not reached. */
  stages: Partial<Record<DmtaStage, string>>;
  /** The furthest stage reached. */
  stage: DmtaStage;
  /** Whole days from designed to tested, when both are dated. */
  cycleDays: number | null;
  sources: string[];
}

const DAY = 24 * 60 * 60 * 1000;

/** Events into one row per compound, in stage order (furthest first), then by compound. A stage is
 * reached when any source says so; its date is the earliest dated event. A later stage implies
 * nothing about an earlier one: a registered compound never designed here has no `designed`. */
export function dmtaRows(events: DmtaEvent[], seriesOf: (compound: string) => string | null = () => null): DmtaRow[] {
  const by = new Map<string, DmtaRow>();
  for (const e of events) {
    let r = by.get(e.compound);
    if (!r)
      by.set(e.compound, r = {compound: e.compound, label: null, series: null, stages: {}, stage: e.stage, cycleDays: null, sources: []});
    if (e.label && !r.label)
      r.label = e.label;
    if (e.series && !r.series)
      r.series = e.series;
    const was = r.stages[e.stage];
    if (was === undefined || (e.on && (was === '' || e.on < was)))
      r.stages[e.stage] = e.on ?? (was ?? '');
    if (!r.sources.includes(e.source))
      r.sources.push(e.source);
  }
  const rows = [...by.values()];
  for (const r of rows) {
    r.series ??= seriesOf(r.compound);
    r.stage = [...DMTA_STAGES].reverse().find((s) => r.stages[s] !== undefined)!;
    const d = r.stages.designed;
    const t = r.stages.tested;
    // whole days: a dataset's date-only value is read at midnight in the reader's zone
    r.cycleDays = d && t ? Math.round((Date.parse(t) - Date.parse(d)) / DAY) : null;
  }
  const order = (s: DmtaStage) => DMTA_STAGES.indexOf(s);
  return rows.sort((a, b) => order(b.stage) - order(a.stage) || a.compound.localeCompare(b.compound));
}

/** Per series: how many compounds reached each stage, and the cycle time designed → tested. */
export interface SeriesCycle {
  series: string;
  compounds: number;
  reached: Record<DmtaStage, number>;
  /** Median days designed → tested over the compounds with both dates; null without any. */
  medianDays: number | null;
  /** How many compounds the median is over. */
  n: number;
}

export function seriesCycles(rows: DmtaRow[]): SeriesCycle[] {
  const groups = new Map<string, DmtaRow[]>();
  for (const r of rows) {
    const k = r.series ?? '—';
    groups.set(k, [...groups.get(k) ?? [], r]);
  }
  return [...groups.entries()].sort(([a], [b]) => a.localeCompare(b)).map(([series, list]) => {
    const days = list.map((r) => r.cycleDays).filter((x): x is number => x !== null).sort((a, b) => a - b);
    const mid = Math.floor(days.length / 2);
    const median = !days.length ? null : days.length % 2 ? days[mid] : (days[mid - 1] + days[mid]) / 2;
    return {series, compounds: list.length, n: days.length, medianDays: median,
      reached: Object.fromEntries(DMTA_STAGES.map((s) => [s, list.filter((r) => r.stages[s] !== undefined).length])) as Record<DmtaStage, number>};
  });
}

/** "KIN-A: 9 compounds · 6 designed, 4 requested, 3 made, 3 tested · cycle 14 days (median of 3)" */
export function describeCycle(c: SeriesCycle): string {
  const reached = DMTA_STAGES.filter((s) => c.reached[s]).map((s) => `${c.reached[s]} ${s}`).join(', ');
  return `${c.series}: ${c.compounds} compound${c.compounds === 1 ? '' : 's'}${reached ? ` · ${reached}` : ''}` +
    (c.medianDays !== null ? ` · cycle ${c.medianDays} days (median of ${c.n})` : ' · no cycle time yet');
}

/** Events a plug-in table delivers: rows with compound, stage and on; a row with an unknown stage
 * or no compound is skipped and counted. */
export function eventsFromRows(rows: Record<string, unknown>[], source: string): {events: DmtaEvent[], skipped: number} {
  const events: DmtaEvent[] = [];
  let skipped = 0;
  for (const r of rows) {
    const compound = r.compound === null || r.compound === undefined ? '' : String(r.compound).trim();
    const stage = String(r.stage ?? '').trim().toLowerCase() as DmtaStage;
    if (!compound || !DMTA_STAGES.includes(stage)) {
      skipped++;
      continue;
    }
    events.push({compound, stage, on: isoOf(r.on), source: r.source ? String(r.source) : source,
      series: r.series ? String(r.series) : null});
  }
  return {events, skipped};
}
