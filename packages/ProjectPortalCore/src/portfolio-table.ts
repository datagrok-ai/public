/** The portfolio as a table: one row per program, with its stage, how long it has been there, its
 * readiness for the next gate, its compounds by role and its work units. Read in a fixed number of
 * batched queries whatever the portfolio's size; only the program dataset of a data condition is
 * read per program, because its SQL selects by `@program`.
 *
 * Three frames, keyed by the program id: the programs, the readiness matrix (programs × the
 * criteria of their next stage) and the compound roles (one row per program and role). */
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, ProgramRow, StageRow} from './generated/db';
import {appliesTo, ContextCriterion, enabledWorkUnits, portfolioByCode, programData, readinessFor, workUnitTable}
  from './portal';
import {CriterionState, nextStage} from './model/gates';
import {DataCondition, DataRow} from './model/conditions';
import {Readiness, ReadinessState} from './model/readiness';
import {COMPOUND_ROLES} from './model/compounds';
import {WORK_UNIT_LABELS, WorkUnitKind} from './model/work-units';
import {existingSpaces} from './workspace';

/** The stage value of a program not placed in its model yet; the board's "Unstaged" column. */
export const UNSTAGED = 'Unstaged';

const LIMIT = 10000;
const BATCH = 100;
const DAY_MS = 24 * 3600 * 1000;

const quoted = (xs: string[]) => xs.map((x) => `"${x.replace(/"/g, '\\"')}"`).join(', ');

async function inBatches<T>(keys: string[], run: (filter: string) => Promise<T[]>): Promise<T[]> {
  const out: T[] = [];
  for (let i = 0; i < keys.length; i += BATCH)
    out.push(...await run(quoted(keys.slice(i, i + BATCH))));
  return out;
}

/** Milliseconds of a Dayjs, a Date or an ISO string; null when absent. */
const ms = (v: any): number | null => {
  if (v === null || v === undefined || v === '')
    return null;
  if (typeof v === 'string')
    return Date.parse(v);
  const n = Number(v.valueOf());
  return Number.isNaN(n) ? null : n;
};

/** The column names of the programs frame, so callers and tests do not spell them twice. */
export const COL = {
  id: 'id', program: 'program', name: 'name', portfolio: 'portfolio', setting: 'setting',
  stage: 'stage', stageName: 'stage name', stageOrdinal: 'stage ordinal', stageInOrder: 'stage in order',
  nextStage: 'next stage',
  status: 'status', modality: 'modality', target: 'target', therapeuticArea: 'therapeutic area',
  source: 'source', sourceStage: 'source stage', enteredStage: 'entered stage', daysInStage: 'days in stage',
  requiredMet: 'required met', requiredTotal: 'required total', readiness: 'readiness', ready: 'ready',
  gaps: 'gaps', compounds: 'compounds', workspace: 'workspace', continued: 'continued',
} as const;

export const roleColumn = (role: string) => `${role} compounds`;
export const unitColumn = (kind: WorkUnitKind) => WORK_UNIT_LABELS[kind].toLowerCase();

/** One program's row before it becomes columns. */
export interface ProgramFacts {
  program: ProgramRow;
  portfolio: PortfolioRow;
  stage: StageRow | null;
  next: StageRow | null;
  enteredOn: number | null;
  daysInStage: number | null;
  readiness: Readiness | null;
  roles: Record<string, number>;
  units: Partial<Record<WorkUnitKind, number>>;
  /** The program's workspace folder exists; a recorded space_id alone does not say so. */
  hasWorkspace?: boolean;
}

export interface PortfolioTables {
  /** One row per program. */
  programs: DG.DataFrame;
  /** One row per program; a column per criterion of a next stage (`<stage>: <criterion>`), its state. */
  readiness: DG.DataFrame;
  /** One row per program and compound role it has. */
  compounds: DG.DataFrame;
  /** The portfolios read, by name. */
  portfolios: PortfolioRow[];
  facts: ProgramFacts[];
}

/** Everything the table needs for [portfolioCode], or for every portfolio the user can see. */
export async function portfolioFacts(portfolioCode?: string, now = Date.now()):
  Promise<{portfolios: PortfolioRow[], stages: Map<string, StageRow[]>, facts: ProgramFacts[]}> {
  const portfolios = portfolioCode ? [await portfolioByCode(portfolioCode)] :
    (await portalDb.portfolios.query({limit: LIMIT})).sort((a, b) => a.name.localeCompare(b.name));
  const portfolioOf = new Map(portfolios.map((p) => [p.id, p]));
  const modelIds = [...new Set(portfolios.map((p) => p.stage_model_id))];
  const [stageRows, programs] = await Promise.all([
    inBatches(modelIds, (b) => portalDb.stages.query({filter: `stage_model_id in (${b})`, limit: LIMIT})),
    inBatches(portfolios.map((p) => p.id), (b) => portalDb.programs.query({filter: `portfolio_id in (${b})`, limit: LIMIT})),
  ]);
  const stages = new Map<string, StageRow[]>(modelIds.map((m) =>
    [m, stageRows.filter((s) => s.stage_model_id === m).sort((a, b) => a.ordinal - b.ordinal)]));
  const stageById = new Map(stageRows.map((s) => [s.id, s]));
  const ids = programs.map((p) => p.id);
  const kinds = [...new Set(portfolios.flatMap(enabledWorkUnits))];

  const [criteria, statuses, transitions, links, ...units] = await Promise.all([
    inBatches(stageRows.map((s) => s.id), (b) => portalDb.gateCriterions.query({filter: `stage_id in (${b})`, limit: LIMIT})),
    inBatches(ids, (b) => portalDb.criterionStatuses.query({filter: `program_id in (${b})`, limit: LIMIT})),
    inBatches(ids, (b) => portalDb.stageTransitions.query({filter: `program_id in (${b})`, limit: LIMIT})),
    inBatches(ids, (b) => portalDb.programCompounds.query({filter: `program_id in (${b})`, limit: LIMIT})),
    ...kinds.map((k) => inBatches(ids, (b) =>
      (workUnitTable(k) as any).query({filter: `program_id in (${b})`, limit: LIMIT}) as Promise<{program_id: string}[]>)),
  ]);

  const group = <T extends {program_id: string}>(rows: T[]) => {
    const m = new Map<string, T[]>();
    for (const r of rows)
      m.set(r.program_id, [...(m.get(r.program_id) ?? []), r]);
    return m;
  };
  const statusesOf = group(statuses);
  const transitionsOf = group(transitions);
  const linksOf = group(links);
  const unitsOf = kinds.map((k, i) => [k, group(units[i] as {program_id: string}[])] as const);

  const facts: ProgramFacts[] = programs.map((program) => {
    const portfolio = portfolioOf.get(program.portfolio_id)!;
    const model = stages.get(program.stage_model_id) ?? stages.get(portfolio.stage_model_id) ?? [];
    const stage = program.stage_id ? stageById.get(program.stage_id) ?? null : null;
    const nextRef = nextStage(model.map((s) => ({id: s.id, code: s.code, name: s.name, ordinal: s.ordinal})),
      program.stage_id ?? null);
    // the stage was entered by the last transition into it; a program never moved entered it when it was made
    const into = (transitionsOf.get(program.id) ?? [])
      .filter((t) => (t.to_stage_id ?? null) === (program.stage_id ?? null) && (t.from_stage_id ?? null) !== (t.to_stage_id ?? null))
      .map((t) => ms(t.decided_on) ?? ms(t.created_on)).filter((t): t is number => t !== null);
    const enteredOn = into.length ? Math.max(...into) : ms(program.created_on);
    const roles: Record<string, number> = {};
    for (const l of linksOf.get(program.id) ?? [])
      roles[l.role] = (roles[l.role] ?? 0) + 1;
    const enabled = enabledWorkUnits(portfolio);
    return {
      program, portfolio, stage, next: nextRef ? stageById.get(nextRef.id) ?? null : null, enteredOn,
      daysInStage: enteredOn === null ? null : Math.max(0, Math.floor((now - enteredOn) / DAY_MS)),
      readiness: null, roles,
      units: Object.fromEntries(unitsOf.filter(([k]) => enabled.includes(k)).map(([k, m]) => [k, m.get(program.id)?.length ?? 0])),
    };
  });

  const live = await existingSpaces(facts.map((f) => f.program.space_id)).catch(() => new Set<string>());
  for (const f of facts)
    f.hasWorkspace = !!f.program.space_id && live.has(f.program.space_id);

  // readiness for the next stage, exactly as the gate panel computes it (readinessFor)
  const criteriaOf = new Map<string, {modality: string | null, criterion: ContextCriterion}[]>();
  for (const c of criteria) {
    criteriaOf.set(c.stage_id, [...(criteriaOf.get(c.stage_id) ?? []), {modality: c.modality ?? null, criterion: {
      id: c.id, stageId: c.stage_id, code: c.code, name: c.name, kind: c.kind, required: c.required !== false,
      rule: c.rule ?? null, condition: (c.condition as DataCondition | undefined) ?? null}}]);
  }
  const work = facts.filter((f) => f.next).map((f) => async () => {
    const mine = (criteriaOf.get(f.next!.id) ?? [])
      .filter((c) => appliesTo(c.modality, f.program.modality)).map((c) => c.criterion);
    const recorded = Object.fromEntries((statusesOf.get(f.program.id) ?? [])
      .map((s) => [s.gate_criterion_id, s.status as CriterionState]));
    let rows: DataRow[] | null = null;
    let dataError: string | undefined;
    if (mine.some((c) => c.condition)) {
      try {
        rows = await programData(f.portfolio, f.program.code);
      } catch (e: any) {
        dataError = e?.message ?? String(e);
      }
    }
    f.readiness = readinessFor({portfolio: f.portfolio, criteria: mine, recorded}, f.next!.id, rows, dataError);
  });
  // the program datasets: a few at a time, not one round trip after another
  for (let i = 0; i < work.length; i += 8)
    await Promise.all(work.slice(i, i + 8).map((w) => w()));

  facts.sort((a, b) => a.portfolio.name.localeCompare(b.portfolio.name) ||
    (a.stage?.ordinal ?? 0) - (b.stage?.ordinal ?? 0) || a.program.code.localeCompare(b.program.code));
  return {portfolios, stages, facts};
}

/** Stage codes in model order, the unstaged first; across models by ordinal, then code. */
function stageOrder(portfolios: PortfolioRow[], stages: Map<string, StageRow[]>): string[] {
  const all = [...new Set(portfolios.map((p) => p.stage_model_id))].flatMap((m) => stages.get(m) ?? []);
  const sorted = portfolios.length === 1 ? all : [...all].sort((a, b) => a.ordinal - b.ordinal || a.code.localeCompare(b.code));
  return [UNSTAGED, ...new Set(sorted.map((s) => s.code))];
}

export const READINESS_COLORS: Record<ReadinessState, string> = {
  met: '#2ca02c', waived: '#9ecae1', in_progress: '#ffbb78', not_met: '#d62728', missing: '#c7c7c7',
};

export function buildTables(portfolios: PortfolioRow[], stages: Map<string, StageRow[]>, facts: ProgramFacts[],
  name = portfolios.length === 1 ? portfolios[0].name : 'All portfolios'): PortfolioTables {
  const col = (type: DG.ColumnType, n: string, f: (x: ProgramFacts) => any) => DG.Column.fromList(type, n, facts.map(f));
  const kinds = [...new Set(portfolios.flatMap(enabledWorkUnits))];
  const S = DG.COLUMN_TYPE.STRING;
  const I = DG.COLUMN_TYPE.INT;
  const programs = DG.DataFrame.fromColumns([
    col(S, COL.id, (x) => x.program.id),
    col(S, COL.program, (x) => x.program.code),
    col(S, COL.name, (x) => x.program.name),
    col(S, COL.portfolio, (x) => x.portfolio.code),
    col(S, COL.setting, (x) => x.portfolio.setting ?? null),
    col(S, COL.stage, (x) => x.stage?.code ?? UNSTAGED),
    col(S, COL.stageName, (x) => x.stage?.name ?? UNSTAGED),
    col(I, COL.stageOrdinal, (x) => x.stage?.ordinal ?? 0),
    // the bar chart sorts its categories by name, not by a column's category order (2026-10-01)
    col(S, COL.stageInOrder, (x) => `${String(x.stage?.ordinal ?? 0).padStart(2, '0')} ${x.stage?.code ?? UNSTAGED}`),
    col(S, COL.nextStage, (x) => x.next?.code ?? null),
    col(S, COL.status, (x) => x.program.status ?? 'active'),
    col(S, COL.modality, (x) => x.program.modality ?? null),
    col(S, COL.target, (x) => x.program.target ?? null),
    col(S, COL.therapeuticArea, (x) => x.program.therapeutic_area ?? null),
    col(S, COL.source, (x) => x.program.source ?? null),
    col(S, COL.sourceStage, (x) => x.program.source_stage ?? null),
    col(DG.COLUMN_TYPE.DATE_TIME, COL.enteredStage, (x) => x.enteredOn === null ? null : new Date(x.enteredOn)),
    col(I, COL.daysInStage, (x) => x.daysInStage),
    col(I, COL.requiredMet, (x) => x.readiness?.requiredMet ?? null),
    col(I, COL.requiredTotal, (x) => x.readiness?.requiredTotal ?? null),
    col(DG.COLUMN_TYPE.FLOAT, COL.readiness, (x) => x.readiness?.score ?? null),
    col(DG.COLUMN_TYPE.BOOL, COL.ready, (x) => x.readiness !== null && x.readiness.score === 1),
    col(S, COL.gaps, (x) => x.readiness?.gaps.join('; ') || null),
    col(I, COL.compounds, (x) => Object.values(x.roles).reduce((a, b) => a + b, 0)),
    ...COMPOUND_ROLES.map((r) => col(I, roleColumn(r), (x) => x.roles[r] ?? 0)),
    ...kinds.map((k) => col(I, unitColumn(k), (x) => x.units[k] ?? null)),
    col(DG.COLUMN_TYPE.BOOL, COL.workspace, (x) => x.hasWorkspace ?? false),
    col(DG.COLUMN_TYPE.BOOL, COL.continued, (x) => !!x.program.predecessor_id),
  ]);
  programs.name = name;
  programs.col(COL.stage)!.setCategoryOrder(stageOrder(portfolios, stages));
  programs.col(COL.readiness)!.meta.format = '0%';

  // the readiness matrix: a column per criterion any program faces at its next stage
  const keys: {key: string, stage: string, code: string}[] = [];
  for (const f of facts) {
    for (const c of f.readiness?.criteria ?? []) {
      const key = `${f.next!.code}: ${c.code}`;
      if (!keys.some((k) => k.key === key))
        keys.push({key, stage: f.next!.code, code: c.code});
    }
  }
  const order = stageOrder(portfolios, stages);
  keys.sort((a, b) => order.indexOf(a.stage) - order.indexOf(b.stage));
  const readiness = DG.DataFrame.fromColumns([
    col(S, COL.id, (x) => x.program.id),
    col(S, COL.program, (x) => x.program.code),
    col(S, COL.portfolio, (x) => x.portfolio.code),
    col(S, COL.nextStage, (x) => x.next?.code ?? null),
    col(DG.COLUMN_TYPE.FLOAT, COL.readiness, (x) => x.readiness?.score ?? null),
    ...keys.map((k) => col(S, k.key, (x) => x.next?.code === k.stage ?
      x.readiness?.criteria.find((c) => c.code === k.code)?.state ?? null : null)),
  ]);
  readiness.name = `${name}: readiness`;
  readiness.col(COL.readiness)!.meta.format = '0%';
  for (const k of keys) {
    const c = readiness.col(k.key)!;
    c.meta.colors.setCategorical(READINESS_COLORS);
  }

  const roleRows = facts.flatMap((f) => COMPOUND_ROLES.filter((r) => f.roles[r]).map((r) => ({f, r, n: f.roles[r]})));
  const compounds = DG.DataFrame.fromColumns([
    DG.Column.fromList(S, COL.id, roleRows.map((x) => x.f.program.id)),
    DG.Column.fromList(S, COL.program, roleRows.map((x) => x.f.program.code)),
    DG.Column.fromList(S, COL.portfolio, roleRows.map((x) => x.f.portfolio.code)),
    DG.Column.fromList(S, COL.stage, roleRows.map((x) => x.f.stage?.code ?? UNSTAGED)),
    DG.Column.fromList(S, 'role', roleRows.map((x) => x.r)),
    DG.Column.fromList(I, COL.compounds, roleRows.map((x) => x.n)),
  ]);
  compounds.name = `${name}: compounds`;
  compounds.col('role')!.setCategoryOrder([...COMPOUND_ROLES]);
  return {programs, readiness, compounds, portfolios, facts};
}

/** The three frames for [portfolioCode], or for every portfolio the user can see. */
export async function portfolioTables(portfolioCode?: string): Promise<PortfolioTables> {
  const {portfolios, stages, facts} = await portfolioFacts(portfolioCode);
  return buildTables(portfolios, stages, facts);
}
