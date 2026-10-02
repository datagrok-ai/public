/** Model recommendation v1 (F-13, FR-51): which catalogued models fit a program, keyed on task,
 * stage, modality and work-unit kind, and which assays the next gate still lacks. A card declares
 * where the model applies (`appliesTo`); a key it leaves out matches anything. Pure: the caller
 * hands in the catalogue, the program's place and its readiness; nothing here touches the platform. */
import {WORK_UNIT_KINDS, WorkUnitKind} from './work-units';
import {CriterionReadiness} from './readiness';
import {columnValue, DataCondition, DataRow, DEFAULT_LEAD_ROLES, holds, isBlank, rowId, rowsInScope} from './conditions';

/** What a model applies to, on its card. Every key is optional; absent or empty means any. */
export interface ModelApplicability {
  /** What the model estimates: potency, ADMET, selectivity, physchem, PK, process, … */
  tasks?: string[];
  /** Stage codes of any stage model: the stages a program is at, or is working towards. */
  stages?: string[];
  /** Program modalities, compared as the criteria compare them (case and punctuation ignored). */
  modalities?: string[];
  /** Work-unit kinds the model works on: series, campaign, study, process_stage. */
  workUnits?: WorkUnitKind[];
  /** Columns of the program dataset the model's outputs estimate: an assay it can stand in for. */
  predicts?: string[];
  /** Order among models that fit equally well; lower first. */
  rank?: number;
}

const KEYS = ['tasks', 'stages', 'modalities', 'workUnits', 'predicts'] as const;

/** Every problem with a card's applicability, as sentences. */
export function validateApplicability(a: ModelApplicability | undefined): string[] {
  if (a === undefined || a === null)
    return [];
  if (typeof a !== 'object' || Array.isArray(a))
    return ['appliesTo must be an object'];
  const errors: string[] = [];
  for (const k of KEYS) {
    const v = (a as any)[k];
    if (v === undefined)
      continue;
    if (!Array.isArray(v) || v.some((x) => typeof x !== 'string' || x.trim() === ''))
      errors.push(`appliesTo.${k} must be a list of non-empty strings`);
  }
  for (const w of a.workUnits ?? []) {
    if (!(WORK_UNIT_KINDS as string[]).includes(w))
      errors.push(`appliesTo.workUnits: "${w}" is not one of ${WORK_UNIT_KINDS.join(', ')}`);
  }
  for (const k of Object.keys(a)) {
    if (!(KEYS as readonly string[]).includes(k) && k !== 'rank')
      errors.push(`appliesTo.${k} is not a key (${[...KEYS, 'rank'].join(', ')})`);
  }
  if (a.rank !== undefined && (typeof a.rank !== 'number' || Number.isNaN(a.rank)))
    errors.push('appliesTo.rank must be a number');
  return errors;
}

// ---------------------------------------------------------------------------------------------
// Missing assays

/** A measurement the gate needs that some rows lack. `column` is null for a criterion that has
 * only a rule in words: the gap is known, the column is not. */
export interface MissingAssay {
  criterion: string;
  criterionName: string;
  required: boolean;
  column: string | null;
  /** Rows that would count once the column is measured. */
  missing: number;
  /** The rows, by compound (or lot, id, code, name). */
  items: string[];
  rule?: string | null;
}

/** Where-columns compared with a number are measurements (`ic50_um < 1`); a column matched against
 * a value or a list (`series_maturity in […]`) classifies a row and is not an assay. */
const measuredBy = (c: DataCondition): string[] => {
  const out = [...(c.measured ?? [])];
  for (const p of c.where ?? []) {
    if (['<', '<=', '>', '>='].includes(p.op) && !out.includes(p.column))
      out.push(p.column);
  }
  return out;
};

/** The assays missing for a stage's open data conditions: per criterion and column, the in-scope
 * rows that lack the column and would otherwise be looked at. A measured column counts the rows
 * passing `where`; a numeric where-column counts the in-scope rows with no value. With [series],
 * only that series' rows. A criterion met or waived needs nothing. */
export function missingAssays(criteria: (CriterionReadiness & {condition?: DataCondition | null, rule?: string | null})[],
  rows: DataRow[] | null, opts: {leadRoles?: string[], series?: string | null} = {}): MissingAssay[] {
  const out: MissingAssay[] = [];
  const inSeries = (row: DataRow) => !opts.series || String(columnValue(row, 'series') ?? '') === opts.series;
  for (const c of criteria) {
    if (c.kind !== 'data_condition' || c.state === 'met' || c.state === 'waived')
      continue;
    const base = {criterion: c.code, criterionName: c.name, required: c.required};
    if (!c.condition || rows === null) {
      if (!c.condition)
        out.push({...base, column: null, missing: 0, items: [], rule: c.rule ?? null});
      continue;
    }
    const scoped = rowsInScope(c.condition.scope ?? 'program', rows, opts.leadRoles ?? DEFAULT_LEAD_ROLES)
      .filter(({row}) => inSeries(row));
    const where = c.condition.where ?? [];
    for (const column of measuredBy(c.condition)) {
      const isMeasured = (c.condition.measured ?? []).includes(column);
      // a measured column: rows that pass where; a where-column: rows that pass the other predicates
      const lacking = scoped.filter(({row}) => isBlank(columnValue(row, column)) &&
        where.filter((p) => isMeasured || p.column !== column).every((p) => holds(row, p)));
      if (lacking.length)
        out.push({...base, column, missing: lacking.length, items: lacking.map(({row, i}) => rowId(row, i))});
    }
  }
  return out;
}

export function describeMissing(m: MissingAssay): string {
  if (m.column === null)
    return `${m.criterionName}: no structured rule${m.rule ? ` (${m.rule})` : ''}; record it by hand`;
  const shown = m.items.slice(0, 5).join(', ') + (m.items.length > 5 ? `, ${m.items.length - 5} more` : '');
  return `${m.column} on ${m.missing} ${m.missing === 1 ? 'row' : 'rows'} (${shown}) for ${m.criterionName}` +
    (m.required ? '' : ' (optional)');
}

// ---------------------------------------------------------------------------------------------
// Recommendation

export interface CandidateModel {
  code: string;
  name: string;
  status?: string | null;
  approval?: string | null;
  /** Portfolios whose catalogue declares the model. */
  portfolios: string[];
  appliesTo?: ModelApplicability | null;
  installed?: boolean;
}

export interface RecommendContext {
  portfolio: string;
  modality?: string | null;
  /** The stage the program is at and the one it works towards (codes); either may be absent. */
  stages: string[];
  /** Work-unit kinds in play: the portfolio's, or the one a series page is about. */
  workUnits: string[];
  /** Narrows to models that declare this task. */
  task?: string | null;
  /** Columns the gate lacks: a model that predicts one ranks first. */
  missingColumns: string[];
}

export interface Recommendation {
  code: string;
  name: string;
  /** Missing columns the model predicts. */
  predicts: string[];
  /** Keys the card declares and the program matches: the more, the more specific the fit. */
  matched: string[];
  reasons: string[];
  approval: string | null;
  installed: boolean;
}

const norm = (s: string | null | undefined) => (s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
const APPROVAL_ORDER: Record<string, number> = {approved: 0, in_review: 1, draft: 2};

/** The models that fit, best first. A candidate is a model the portfolio's catalogue declares, or
 * any model whose card says where it applies; a deprecated model or a retired card never is. Each
 * key the card declares must match: stage (current or target), modality, work-unit kind, and the
 * task when one is asked for. Then: models predicting a missing column first, then the more
 * specific fit, approved before in review before draft, the card's rank, the name. */
export function recommend(models: CandidateModel[], ctx: RecommendContext): Recommendation[] {
  const out: (Recommendation & {rank: number})[] = [];
  for (const m of models) {
    const a = m.appliesTo ?? {};
    const declared = (k: keyof ModelApplicability) => Array.isArray(a[k]) && (a[k] as string[]).length > 0;
    if (m.status === 'deprecated' || m.approval === 'retired')
      continue;
    if (!m.portfolios.includes(ctx.portfolio) && !KEYS.some(declared))
      continue;
    const matched: string[] = [];
    const reasons: string[] = [];
    if (declared('stages')) {
      const hit = ctx.stages.find((s) => a.stages!.some((x) => norm(x) === norm(s)));
      if (!hit)
        continue;
      matched.push('stage');
      reasons.push(`stage ${hit}`);
    }
    if (declared('modalities')) {
      if (!ctx.modality || !a.modalities!.some((x) => norm(x) === norm(ctx.modality)))
        continue;
      matched.push('modality');
      reasons.push(ctx.modality);
    }
    if (declared('workUnits')) {
      const hit = ctx.workUnits.find((w) => (a.workUnits as string[]).includes(w));
      if (!hit)
        continue;
      matched.push('work unit');
      reasons.push(hit.replace('_', ' '));
    }
    if (ctx.task) {
      if (!declared('tasks') || !a.tasks!.some((t) => norm(t) === norm(ctx.task)))
        continue;
      matched.push('task');
    }
    if (declared('tasks'))
      reasons.push(a.tasks!.join(', '));
    const predicts = (a.predicts ?? []).filter((c) => ctx.missingColumns.some((x) => norm(x) === norm(c)));
    if (predicts.length)
      reasons.push(`predicts ${predicts.join(', ')}, missing for the gate`);
    if (m.portfolios.includes(ctx.portfolio) && !matched.length)
      reasons.push(`in the ${ctx.portfolio} catalogue`);
    out.push({code: m.code, name: m.name, predicts, matched, reasons, approval: m.approval ?? null,
      installed: m.installed !== false, rank: a.rank ?? Number.MAX_SAFE_INTEGER});
  }
  out.sort((x, y) => y.predicts.length - x.predicts.length || y.matched.length - x.matched.length ||
    (APPROVAL_ORDER[x.approval ?? ''] ?? 3) - (APPROVAL_ORDER[y.approval ?? ''] ?? 3) || x.rank - y.rank ||
    x.name.localeCompare(y.name));
  return out.map(({rank: _rank, ...r}) => r);
}

/** For each missing column, the recommended models that predict it. */
export function predictorsOf(missing: MissingAssay[], recs: Recommendation[]): Map<string, Recommendation[]> {
  const out = new Map<string, Recommendation[]>();
  for (const m of missing) {
    if (m.column && !out.has(m.column))
      out.set(m.column, recs.filter((r) => r.predicts.some((p) => norm(p) === norm(m.column))));
  }
  return out;
}

/** Values a model's dialog opens with: a dataframe input gets the program's data, an input named
 * program, portfolio or series gets that code. Anything else keeps the function's default. */
export function prefillInputs(inputs: {name: string, type: string}[],
  ctx: {portfolio: string, program: string, series?: string | null}): {values: Record<string, unknown>, needsData: string[]} {
  const values: Record<string, unknown> = {};
  const needsData: string[] = [];
  for (const p of inputs) {
    const n = p.name.toLowerCase();
    if (p.type === 'dataframe')
      needsData.push(p.name);
    else if (p.type === 'string' && n === 'program')
      values[p.name] = ctx.program;
    else if (p.type === 'string' && n === 'portfolio')
      values[p.name] = ctx.portfolio;
    else if (p.type === 'string' && n === 'series' && ctx.series)
      values[p.name] = ctx.series;
  }
  return {values, needsData};
}
