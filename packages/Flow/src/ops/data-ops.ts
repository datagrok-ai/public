/** Flow-native row/column/aggregation verbs: the core originals take predicate/spec types no canvas
 *  socket can feed; these take primitives, delegate to the platform engines, and return new tables. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Scratch name for the condition column — only ever exists detached (`addColumn: false`). */
const CONDITION_COLUMN = '~flow condition';

/** Literal strings, not `DG.AGG` references — the annotation generator only reads literal arrays. */
export const AGGREGATION_TYPES = [
  'count', 'values', 'unique', 'nulls',
  'min', 'max', 'sum', 'avg', 'med', 'stdev', 'variance', 'q1', 'q2', 'q3', 'skew', 'kurt',
  'first', 'concat unique',
];

/** `count` is the only aggregation meaningful with no column at all. */
const COLUMNLESS_AGGREGATIONS = new Set(['count']);

/** One entry of the `aggregations` JSON. `column` may be blank only for COLUMNLESS_AGGREGATIONS. */
export interface AggregationSpec {
  column: string;
  type: string;
  name?: string;
}

/** Rows matching `condition`, via `core:AddNewColumn` with aux `addColumn: false` —
 *  computes the formula column WITHOUT appending it, so the user's table is never touched. */
export async function conditionMask(table: DG.DataFrame, condition: string): Promise<DG.BitSet> {
  const expression = String(condition ?? '').trim();
  if (expression === '')
    throw new Error('Condition is empty — write a boolean expression, e.g. ${age} > 30');

  // `auto`, deliberately NOT `bool`: asking for a bool COERCES — `${score} + 1` comes back
  // as a bool column of nulls and silently matches nothing.
  const call = DG.Func.byName('AddNewColumn').prepare({
    table, expression, name: CONDITION_COLUMN, type: 'auto',
    treatAsString: false, subscribeOnChanges: false,
  });
  call.setAuxValue('addColumn', false);
  await call.call();

  const column = call.getOutputParamValue() as DG.Column | null;
  if (!column)
    throw new Error(`Condition "${expression}" produced no values — check the column names`);
  if (column.type !== DG.COLUMN_TYPE.BOOL) {
    throw new Error(`Condition "${expression}" is ${column.type}, not a true/false expression. ` +
      'Compare something, e.g. ${age} > 30');
  }
  // A null (a row the formula could not evaluate) must not silently pass a filter.
  return DG.BitSet.create(table.rowCount, (i) => column.get(i) === true);
}

export async function filterRows(table: DG.DataFrame, condition: string): Promise<DG.DataFrame> {
  return named(table.clone(await conditionMask(table, condition)), `${table.name} filtered`);
}

/** The complement of {@link filterRows}. */
export async function deleteRows(table: DG.DataFrame, condition: string): Promise<DG.DataFrame> {
  const mask = await conditionMask(table, condition);
  return named(table.clone(mask.invert()), `${table.name} trimmed`);
}

/** Rows matching the condition, keeping only the chosen columns; none chosen keeps all. */
export async function extractRows(
  table: DG.DataFrame, condition: string, columns?: DG.Column[] | null,
): Promise<DG.DataFrame> {
  const names = columnNames(columns);
  const mask = await conditionMask(table, condition);
  return named(table.clone(mask, names.length > 0 ? names : null), `${table.name} extract`);
}

/** Deliberately a mutation — selection is table state; a copy would carry it where nobody sees. */
export async function selectRows(
  table: DG.DataFrame, condition: string, clearSelection: boolean = true,
): Promise<DG.DataFrame> {
  const mask = await conditionMask(table, condition);
  if (clearSelection) table.selection.setAll(false, false);
  table.selection.or(mask);
  return table;
}

/** Drawn without replacement from a seeded generator so a pipeline replays identically. Pure. */
export function randomIndices(rowCount: number, count: number, seed: number): number[] {
  const n = Math.max(0, Math.min(Math.floor(count), rowCount));
  const pool = Array.from({length: rowCount}, (_, i) => i);
  const random = seededRandom(seed);
  // Partial Fisher-Yates: only the first `n` slots need to be settled.
  for (let i = 0; i < n; i++) {
    const j = i + Math.floor(random() * (rowCount - i));
    const t = pool[i];
    pool[i] = pool[j];
    pool[j] = t;
  }
  return pool.slice(0, n).sort((a, b) => a - b);
}

export function filterRandomRows(table: DG.DataFrame, count: number, seed: number = 42): DG.DataFrame {
  const picked = new Set(randomIndices(table.rowCount, count, seed));
  return named(table.clone(DG.BitSet.create(table.rowCount, (i) => picked.has(i))),
    `${table.name} sample`);
}

export function selectRandomRows(
  table: DG.DataFrame, count: number, seed: number = 42, clearSelection: boolean = true,
): DG.DataFrame {
  const picked = new Set(randomIndices(table.rowCount, count, seed));
  if (clearSelection) table.selection.setAll(false, false);
  for (const i of picked) table.selection.set(i, true, false);
  table.selection.fireChanged();
  return table;
}

/** A copy of the table without the chosen columns. */
export function deleteColumns(table: DG.DataFrame, columns: DG.Column[]): DG.DataFrame {
  const dropped = new Set(columnNames(columns));
  if (dropped.size === 0)
    throw new Error('Pick at least one column to delete');
  const kept = table.columns.names().filter((n) => !dropped.has(n));
  if (kept.length === 0)
    throw new Error('Deleting every column would leave nothing — keep at least one');
  return named(table.clone(null, kept), table.name);
}

/** Tags are column metadata, so there is nothing to copy — the mutation IS the result. */
export function tagColumns(
  table: DG.DataFrame, columns: DG.Column[], tag: string, value: string,
): DG.DataFrame {
  const key = String(tag ?? '').trim();
  if (key === '')
    throw new Error('Tag name is empty');
  const names = columnNames(columns);
  if (names.length === 0)
    throw new Error('Pick at least one column to tag');
  for (const name of names) {
    const column = table.col(name);
    if (!column)
      throw new Error(`Column "${name}" is not in ${table.name}`);
    column.setTag(key, String(value ?? ''));
  }
  return table;
}

/** `core:AddNewColumn` with Flow's inline formula editor and none of the reactivity plumbing. */
export async function expressionToColumn(
  table: DG.DataFrame, expression: string, name: string, type: string = 'auto',
): Promise<DG.Column> {
  const formula = String(expression ?? '').trim();
  if (formula === '')
    throw new Error('Expression is empty');
  const columnName = String(name ?? '').trim() || table.columns.getUnusedName('result');
  const column = await grok.functions.call('AddNewColumn', {
    table, expression: formula, name: columnName, type,
  }) as DG.Column | null;
  if (!column)
    throw new Error(`Expression "${formula}" produced no column`);
  return column;
}

/** Tolerates a blank or half-typed value (runs while the panel renders); unknown
 *  aggregation types are dropped — {@link aggregationProblems} reports them. */
export function parseAggregations(stored: unknown): AggregationSpec[] {
  if (stored === null || stored === undefined || String(stored).trim() === '') return [];
  let parsed: unknown;
  try {
    parsed = typeof stored === 'string' ? JSON.parse(String(stored)) : stored;
  } catch {
    return [];
  }
  if (!Array.isArray(parsed)) return [];
  return parsed
    .filter((e): e is Record<string, unknown> => !!e && typeof e === 'object')
    .map((e) => ({
      column: String(e.column ?? ''),
      type: String(e.type ?? ''),
      name: e.name === undefined || e.name === null ? undefined : String(e.name),
    }))
    .filter((a) => AGGREGATION_TYPES.includes(a.type));
}

/** Shared by the node's readiness check and the function, so hint and runtime failure agree. */
export function aggregationProblems(specs: readonly AggregationSpec[]): string[] {
  if (specs.length === 0) return ['at least one aggregation'];
  const missing = specs
    .filter((a) => !COLUMNLESS_AGGREGATIONS.has(a.type) && a.column.trim() === '')
    .map((a) => `a column for "${a.type}"`);
  return Array.from(new Set(missing));
}

/** A pivot column turns the group-by into a pivot table — hence also findable as "pivot". */
export function aggregate(
  table: DG.DataFrame,
  groupByColumns: DG.Column[] | null,
  aggregations: string,
  pivotColumns?: DG.Column[] | null,
): DG.DataFrame {
  const specs = parseAggregations(aggregations);
  const problems = aggregationProblems(specs);
  if (problems.length > 0)
    throw new Error(`Aggregate needs ${problems.join(', ')}`);

  const builder = table.groupBy(columnNames(groupByColumns));
  for (const spec of specs) {
    const column = COLUMNLESS_AGGREGATIONS.has(spec.type) ? null : spec.column;
    // The builder's TS signature is the `AGG` union, but Dart takes the name as a
    // string — how `concat unique` (a `STR_AGG`, not an `AGG`) reaches it.
    builder.add(spec.type as DG.AggregationType, column, spec.name?.trim() || null);
  }
  for (const name of columnNames(pivotColumns))
    builder.pivot(name);
  return named(builder.aggregate(), `${table.name} aggregated`);
}

/** Wide → long; delegates to `core:Unpivot`, whose only canvas flaw is bare string-list params. */
export async function unpivot(
  table: DG.DataFrame,
  copyColumns: DG.Column[] | null,
  mergeColumns: DG.Column[],
  categoryColumnName: string = 'Category',
  valueColumnName: string = 'Value',
): Promise<DG.DataFrame> {
  const merge = columnNames(mergeColumns);
  if (merge.length === 0)
    throw new Error('Pick at least one column to merge — those become the (category, value) rows');
  const result = await grok.functions.call('Unpivot', {
    table,
    copyColumns: columnNames(copyColumns),
    mergeColumns: merge,
    categoryColumnName: String(categoryColumnName ?? '').trim() || 'Category',
    valueColumnName: String(valueColumnName ?? '').trim() || 'Value',
  }) as DG.DataFrame;
  return named(result, `${table.name} unpivoted`);
}

/** Tolerates the null/undefined a blank `column_list` compiles to. */
function columnNames(columns: DG.Column[] | null | undefined): string[] {
  if (!columns) return [];
  return columns.filter((c) => c != null).map((c) => c.name);
}

/** So output tabs and previews read as something other than `table`. */
function named(table: DG.DataFrame, name: string): DG.DataFrame {
  table.name = name;
  return table;
}

/** mulberry32 — `Math.random` would make re-runs non-reproducible and break invalidation. */
function seededRandom(seed: number): () => number {
  let a = (Math.floor(seed) || 0) >>> 0;
  return () => {
    a = (a + 0x6D2B79F5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
