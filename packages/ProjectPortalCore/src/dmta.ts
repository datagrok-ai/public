/** DMTA status (F-24): reads each stage's source into events and hands them to the pure rules
 * (src/model/dmta.ts). Sources are pluggable: a package may register one in code
 * (`registerDmtaSource`), and a portfolio may name functions in `dmta.sources`, until request
 * objects (F-49) exist and become a source of their own. */
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {compoundsOf, programData} from './portal';
import {designSetsOf} from './design-sets';
import {observationsOf} from './observations';
import {datasetsOf, programTable} from './program-data';
import {MAIN_DATASET} from './model/program-data';
import {DmtaConfig, DmtaEvent, DmtaRow, dmtaRows, DMTA_STAGES, eventsFromRows, proposalCompound, SeriesCycle, seriesCycles,
  testedColumns, testedEvents} from './model/dmta';
import {filterValue, ROW_LIMIT} from './rows';

const iso = (d: any): string | null => d ? dayjs(d).toISOString() : null;

export interface DmtaContext {
  portfolio: PortfolioRow;
  program: ProgramRow;
}

/** A source of events, beside the core's own. */
export interface DmtaSource {
  name: string;
  events(ctx: DmtaContext): Promise<DmtaEvent[]>;
}

const registered: DmtaSource[] = [];

/** Adds a source of DMTA events for every program, e.g. a request system's status (F-49). A source
 * of the same name is replaced. */
export function registerDmtaSource(source: DmtaSource): void {
  const i = registered.findIndex((s) => s.name === source.name);
  if (i >= 0)
    registered.splice(i, 1, source);
  else
    registered.push(source);
}

export function unregisterDmtaSource(name: string): void {
  const i = registered.findIndex((s) => s.name === name);
  if (i >= 0)
    registered.splice(i, 1);
}

/** Designed and requested: each proposal of each design set; it is the compound it became, when
 * registered, else its label. */
async function designEvents(ctx: DmtaContext, compounds: {code: string, smiles?: string | null}[]): Promise<DmtaEvent[]> {
  const out: DmtaEvent[] = [];
  for (const {set, proposals} of await designSetsOf(ctx.program.id)) {
    const source = `design set ${set.code}`;
    for (const p of proposals) {
      const compound = proposalCompound(p, compounds) ?? p.label;
      out.push({compound, stage: 'designed', on: iso(set.designed_on ?? set.created_on), source, series: set.series_code, label: p.label});
      if (p.chosen)
        out.push({compound, stage: 'requested', on: iso(p.chosen_on), source, series: set.series_code, label: p.label});
    }
  }
  return out;
}

/** Tested: the first dated row per compound of the configured dataset, else of the first dataset
 * with a compound and a date column (further datasets before the program dataset). */
async function testedEventsOf(ctx: DmtaContext, config: DmtaConfig | null, notes: string[]): Promise<DmtaEvent[]> {
  const all = datasetsOf(ctx.portfolio);
  const tried = config?.tested ? all.filter((d) => d.code === config.tested!.dataset) :
    [...all.filter((d) => d.code !== MAIN_DATASET), ...all.filter((d) => d.code === MAIN_DATASET)];
  for (const d of tried) {
    try {
      const t = await programTable(ctx.portfolio.code, ctx.program.code, d.code);
      const cols = testedColumns(t.df.columns.names(), config?.tested);
      if (!cols)
        continue;
      const rows = Array.from({length: t.df.rowCount}, (_, i) => ({[cols.compound]: t.df.get(cols.compound, i),
        [cols.date]: t.df.get(cols.date, i)}));
      return testedEvents(rows, cols, d.name);
    } catch (e: any) {
      notes.push(`${d.name} could not be read: ${e?.message ?? e}`);
    }
  }
  notes.push(config?.tested ? `${config.tested.dataset} has no compound and date column` :
    'No dataset has a compound and a date column: nothing is marked tested');
  return [];
}

/** Analyzed: a current observation on the compound, or its role set by hand. */
async function analyzedEvents(ctx: DmtaContext, compoundIds: Map<string, string>): Promise<DmtaEvent[]> {
  const out: DmtaEvent[] = (await observationsOf(ctx.program.id).catch(() => []))
    .filter((o) => o.subject === 'compound' && o.code && o.status === 'current')
    .map((o) => ({compound: o.code!, stage: 'analyzed' as const, on: o.on, source: 'observation'}));
  const changes = await portalDb.compoundRoleChanges.query({filter: `program_id = ${filterValue(ctx.program.id)}`, limit: ROW_LIMIT})
    .catch(() => []);
  for (const c of changes.filter((x) => x.source === 'manual')) {
    const code = compoundIds.get(c.compound_id);
    if (code)
      out.push({compound: code, stage: 'analyzed', on: iso(c.changed_on ?? c.created_on), source: 'role set by hand'});
  }
  return out;
}

/** Events from the portfolio's `dmta.sources` functions. */
async function pluginEvents(ctx: DmtaContext, config: DmtaConfig | null, notes: string[]): Promise<DmtaEvent[]> {
  const out: DmtaEvent[] = [];
  for (const name of config?.sources ?? []) {
    try {
      const [pkg, fn] = name.split(':');
      const func = DG.Func.find({package: pkg, name: fn})[0];
      if (!func) {
        notes.push(`${name} is not installed here`);
        continue;
      }
      const params: Record<string, unknown> = {};
      for (const p of func.inputs) {
        if (p.name === 'portfolio')
          params.portfolio = ctx.portfolio.code;
        if (p.name === 'program')
          params.program = ctx.program.code;
      }
      const df = await func.apply(params);
      if (!(df instanceof DG.DataFrame)) {
        notes.push(`${name} returned no table`);
        continue;
      }
      const rows = Array.from({length: df.rowCount}, (_, i) => Object.fromEntries(df.columns.names().map((c) => [c, df.get(c, i)])));
      const r = eventsFromRows(rows, name);
      out.push(...r.events);
      if (r.skipped)
        notes.push(`${name}: ${r.skipped} row(s) without a compound or a known stage (${DMTA_STAGES.join(', ')})`);
    } catch (e: any) {
      notes.push(`${name} failed: ${e?.message ?? e}`);
    }
  }
  for (const s of registered) {
    try {
      out.push(...await s.events(ctx));
    } catch (e: any) {
      notes.push(`${s.name} failed: ${e?.message ?? e}`);
    }
  }
  return out;
}

export interface DmtaStatus {
  rows: DmtaRow[];
  cycles: SeriesCycle[];
  /** What could not be read, or was skipped. */
  notes: string[];
}

/** Where each compound of the program stands, and the cycle time per series. */
export async function dmtaOf(portfolio: PortfolioRow, program: ProgramRow): Promise<DmtaStatus> {
  const ctx = {portfolio, program};
  const config = (portfolio.dmta ?? null) as DmtaConfig | null;
  const notes: string[] = [];
  const linked = await compoundsOf(program.id).catch(() => []);
  const compounds = linked.map((c) => ({code: String(c.compound.code), smiles: c.compound.smiles ?? null}));
  // made: registered in the program (the registry projection), date unknown
  const made: DmtaEvent[] = linked.filter((c) => !c.link.missing_from_feed)
    .map((c) => ({compound: String(c.compound.code), stage: 'made' as const, on: null, source: 'registry'}));
  const events = [
    ...await designEvents(ctx, compounds),
    ...made,
    ...await testedEventsOf(ctx, config, notes),
    ...await analyzedEvents(ctx, new Map(linked.map((c) => [c.link.compound_id, String(c.compound.code)]))),
    ...await pluginEvents(ctx, config, notes),
  ];
  // a compound's series from the program dataset when no design set says it
  const data = await programData(portfolio, program.code).catch(() => null);
  const series = new Map((data ?? []).filter((r) => r.compound && r.series).map((r) => [String(r.compound), String(r.series)]));
  const rows = dmtaRows(events, (c) => series.get(c) ?? null);
  return {rows, cycles: seriesCycles(rows), notes};
}

const day = (s: string | undefined) => s === undefined ? '' : s === '' ? 'yes' : s.slice(0, 10);

export function dmtaFrame(s: DmtaStatus, name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Compound', s.rows.map((r) => r.compound)),
    DG.Column.fromStrings('Design label', s.rows.map((r) => r.label ?? '')),
    DG.Column.fromStrings('Series', s.rows.map((r) => r.series ?? '')),
    DG.Column.fromStrings('Stage', s.rows.map((r) => r.stage)),
    ...DMTA_STAGES.map((st) => DG.Column.fromStrings(st[0].toUpperCase() + st.slice(1), s.rows.map((r) => day(r.stages[st])))),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Cycle days', s.rows.map((r) => r.cycleDays)),
    DG.Column.fromStrings('Sources', s.rows.map((r) => r.sources.join('; '))),
  ]);
  df.col('Stage')!.setCategoryOrder(DMTA_STAGES);
  df.name = name;
  return df;
}

export function cyclesFrame(s: DmtaStatus, name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Series', s.cycles.map((c) => c.series)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Compounds', s.cycles.map((c) => c.compounds)),
    ...DMTA_STAGES.map((st) => DG.Column.fromList(DG.COLUMN_TYPE.INT, st[0].toUpperCase() + st.slice(1), s.cycles.map((c) => c.reached[st]))),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Median cycle days', s.cycles.map((c) => c.medianDays)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Cycles measured', s.cycles.map((c) => c.n)),
  ]);
  df.name = name;
  return df;
}
