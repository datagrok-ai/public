/** Data conditions: the one rule language gate criteria use over a program's data (07 §2.10).
 * A condition is a declarative predicate over the rows of the program dataset (one row per
 * compound, lot or other unit the portfolio's query returns): which rows are in scope, which of
 * them qualify, which columns must be measured, and how many must qualify. Pure: the caller runs
 * the program dataset query and hands the rows in, so the same rule is checked in a test, in the
 * gate panel and before a decision is recorded. */
import type {DatasetViewer} from './program-data';

export type ConditionScope = 'program' | 'lead_set' | 'series';
export type ConditionOp = '<' | '<=' | '>' | '>=' | '=' | '!=' | 'in' | 'not_in' | 'present' | 'absent';

export const CONDITION_OPS: ConditionOp[] = ['<', '<=', '>', '>=', '=', '!=', 'in', 'not_in', 'present', 'absent'];

export interface Predicate {
  column: string;
  op: ConditionOp;
  value?: number | string | boolean | (number | string)[];
}

export interface DataCondition {
  /** Which rows count: every row (`program`), rows whose role is a lead role (`lead_set`), or rows
   * grouped by their `series` column (`series`). Default `program`. */
  scope?: ConditionScope;
  /** Rows qualify when every predicate holds. */
  where?: Predicate[];
  /** Columns a qualifying row must have a value in: an assay panel. */
  measured?: string[];
  /** How many qualifying rows are needed; per series when the scope is `series`. Default 1. */
  minCount?: number;
  /** With scope `series`: how many series must reach `minCount`. Default 1. */
  minSeries?: number;
}

/** The program dataset query of a portfolio: one SQL on a named connection, with `@program`
 * standing for the program's code. */
export interface ProgramDataConfig {
  connection: string;
  sql: string;
  /** Values of the `role` column that make a row part of the lead set. */
  leadRoles?: string[];
  /** What the program page calls it ("SAR table"); default "Program data". */
  name?: string;
  /** Viewers to open it with; default: by the portfolio's setting and the columns. */
  viewers?: DatasetViewer[];
}

export const DEFAULT_LEAD_ROLES = ['lead', 'backup', 'candidate'];

export type DataRow = Record<string, unknown>;

export interface ConditionResult {
  met: boolean;
  /** Rows the scope admits. */
  inScope: number;
  /** Rows that pass `where` and have every `measured` column. */
  qualifying: number;
  required: number;
  /** With scope `series`: qualifying rows per series, and the series that reach `minCount`. */
  series?: Record<string, number>;
  seriesMet?: string[];
  /** Per `measured` column: rows passing `where` that have it. */
  coverage: Record<string, number>;
  /** Rows passing `where` that miss at least one `measured` column, with what they miss. */
  gaps: {id: string, missing: string[]}[];
  /** One sentence for the gate panel and the decision log. */
  summary: string;
}

const isBlank = (v: unknown) => v === null || v === undefined || v === '' ||
  (typeof v === 'number' && Number.isNaN(v));

const num = (v: unknown): number | null => {
  if (typeof v === 'number')
    return Number.isNaN(v) ? null : v;
  if (typeof v === 'string' && v.trim() !== '' && !Number.isNaN(Number(v)))
    return Number(v);
  return null;
};

const same = (a: unknown, b: unknown): boolean => {
  if (isBlank(a) || isBlank(b))
    return false;
  const x = num(a);
  const y = num(b);
  if (x !== null && y !== null)
    return x === y;
  return String(a).toLowerCase() === String(b).toLowerCase();
};

/** Column names are matched case-insensitively: query engines differ in what case they return. */
function get(row: DataRow, column: string): unknown {
  if (column in row)
    return row[column];
  const key = Object.keys(row).find((k) => k.toLowerCase() === column.toLowerCase());
  return key === undefined ? undefined : row[key];
}

export function holds(row: DataRow, p: Predicate): boolean {
  const v = get(row, p.column);
  switch (p.op) {
  case 'present': return !isBlank(v);
  case 'absent': return isBlank(v);
  case 'in':
  case 'not_in': {
    const list = Array.isArray(p.value) ? p.value : [p.value];
    const hit = list.some((x) => same(v, x));
    return p.op === 'in' ? hit : !isBlank(v) && !hit;
  }
  case '=': return same(v, p.value);
  case '!=': return !isBlank(v) && !same(v, p.value);
  default: {
    const x = num(v);
    const y = num(p.value);
    if (x === null || y === null)
      return false;
    return p.op === '<' ? x < y : p.op === '<=' ? x <= y : p.op === '>' ? x > y : x >= y;
  }
  }
}

/** Every problem with a condition, as sentences; empty when it can be evaluated. */
export function validateCondition(c: DataCondition): string[] {
  const errors: string[] = [];
  if (c.scope !== undefined && !['program', 'lead_set', 'series'].includes(c.scope))
    errors.push(`scope "${c.scope}" is not one of program, lead_set, series`);
  for (const p of c.where ?? []) {
    if (!p.column)
      errors.push('a predicate has no column');
    if (!(CONDITION_OPS as string[]).includes(p.op))
      errors.push(`predicate on ${p.column}: operator "${p.op}" is not one of ${CONDITION_OPS.join(' ')}`);
    else if (['<', '<=', '>', '>='].includes(p.op) && num(p.value) === null)
      errors.push(`predicate on ${p.column}: "${p.op}" needs a number`);
    else if (['in', 'not_in'].includes(p.op) && !Array.isArray(p.value))
      errors.push(`predicate on ${p.column}: "${p.op}" needs a list`);
    else if (['=', '!='].includes(p.op) && (p.value === undefined || p.value === null))
      errors.push(`predicate on ${p.column}: "${p.op}" needs a value`);
  }
  for (const m of c.measured ?? []) {
    if (!m)
      errors.push('a measured column has no name');
  }
  if (c.minCount !== undefined && (!Number.isInteger(c.minCount) || c.minCount < 1))
    errors.push(`minCount ${c.minCount} is not a positive integer`);
  if (c.minSeries !== undefined) {
    if (c.scope !== 'series')
      errors.push('minSeries applies only to scope "series"');
    else if (!Number.isInteger(c.minSeries) || c.minSeries < 1)
      errors.push(`minSeries ${c.minSeries} is not a positive integer`);
  }
  if (!c.where?.length && !c.measured?.length && (c.minCount ?? 1) <= 1 && c.scope !== 'series')
    errors.push('the condition holds for any program with one row of data; give it where, measured or minCount');
  return errors;
}

const rowId = (row: DataRow, i: number): string => {
  for (const k of ['compound', 'id', 'code', 'lot', 'name']) {
    const v = get(row, k);
    if (!isBlank(v))
      return String(v);
  }
  return `row ${i + 1}`;
};

const plural = (n: number, one: string, many = `${one}s`) => `${n} ${n === 1 ? one : many}`;

/** Evaluates [c] over [rows], the program's dataset. */
export function evaluateCondition(c: DataCondition, rows: DataRow[], leadRoles = DEFAULT_LEAD_ROLES): ConditionResult {
  const scope = c.scope ?? 'program';
  const required = c.minCount ?? 1;
  const roles = leadRoles.map((r) => r.toLowerCase());
  const indexed = rows.map((row, i) => ({row, i}));
  const scoped = scope === 'lead_set' ?
    indexed.filter(({row}) => roles.includes(String(get(row, 'role') ?? '').toLowerCase())) :
    scope === 'series' ? indexed.filter(({row}) => !isBlank(get(row, 'series'))) : indexed;
  const passing = scoped.filter(({row}) => (c.where ?? []).every((p) => holds(row, p)));
  const measured = c.measured ?? [];
  const coverage: Record<string, number> = {};
  for (const m of measured)
    coverage[m] = passing.filter(({row}) => !isBlank(get(row, m))).length;
  const gaps: {id: string, missing: string[]}[] = [];
  const qualifying = passing.filter(({row, i}) => {
    const missing = measured.filter((m) => isBlank(get(row, m)));
    if (missing.length > 0)
      gaps.push({id: rowId(row, i), missing});
    return missing.length === 0;
  });

  const scopeWord = scope === 'lead_set' ? 'lead-set rows' : scope === 'series' ? 'rows in a series' : 'rows';
  const what = measured.length ? ` with ${measured.join(', ')}` : '';
  if (scope === 'series') {
    const series: Record<string, number> = {};
    for (const {row} of qualifying) {
      const s = String(get(row, 'series'));
      series[s] = (series[s] ?? 0) + 1;
    }
    const minSeries = c.minSeries ?? 1;
    const seriesMet = Object.entries(series).filter(([, n]) => n >= required).map(([s]) => s).sort();
    return {
      met: seriesMet.length >= minSeries, inScope: scoped.length, qualifying: qualifying.length, required,
      series, seriesMet, coverage, gaps,
      summary: `${plural(seriesMet.length, 'series', 'series')} of ${minSeries} needed have ≥${required} ` +
        `qualifying ${required === 1 ? 'row' : 'rows'}${what}` +
        (seriesMet.length ? ` (${seriesMet.join(', ')})` : ''),
    };
  }
  return {
    met: qualifying.length >= required, inScope: scoped.length, qualifying: qualifying.length, required,
    coverage, gaps,
    summary: `${qualifying.length} of ${required} needed${what}; ${plural(scoped.length, scopeWord.replace(/s$/, ''), scopeWord)} in scope`,
  };
}
