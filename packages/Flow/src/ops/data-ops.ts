/** Flow-native data operations — the row/column/aggregation verbs a pipeline
 *  needs, with signatures a canvas can actually wire.
 *
 *  **Why these exist.** The platform already has `core:FilterRows`,
 *  `core:SelectRows`, `core:DeleteRows`, `core:ExtractRows`, `core:Subset`,
 *  `core:DeleteColumns`, `core:TagColumns`, `core:Aggregate` and friends — but
 *  every one of them takes a `TableRowFilterCall` / `ColFilterCall` (a nested
 *  *function call* describing the predicate; see
 *  `ddt/lib/src/functions/rows_predicates.dart`). On a canvas that type has no
 *  editor, no socket anything can feed, and no meaningful default, so the whole
 *  family had to be left out of the catalog. `core:Aggregate` is the same story
 *  one level up: `List<GroupAggregation>` / `List<TableJoin>` /
 *  `List<FieldPredicate>` are JSON-ish object lists, and `core:Unpivot` takes
 *  bare `string_list`s where a column picker belongs.
 *
 *  These take primitives instead — a table, real `column`/`column_list` slots
 *  (which get Flow's column picker for free), constrained `choices`, and a
 *  condition as a plain string edited with the platform's own formula editor
 *  (see `panel/editors/expression-editor.ts`).
 *
 *  **Nothing here reimplements an engine.** Conditions are evaluated by
 *  `core:AddNewColumn` — the same formula engine the editor validates against,
 *  so anything the editor accepts runs. Aggregation is `DataFrame.groupBy`.
 *  Unpivot delegates to `core:Unpivot`.
 *
 *  **Everything returns a new table** rather than mutating in place. Core's
 *  originals are in-place because they drive a grid; a pipeline step whose
 *  result is "the input, but different now" gives the canvas nothing to wire
 *  onward and makes re-running non-idempotent. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Scratch name for the condition column — never appended to the user's table
 *  (`addColumn: false` below), so it only ever exists as a detached column. */
const CONDITION_COLUMN = '~flow condition';

/** Aggregations offered by {@link aggregate}, in the order the editor lists
 *  them. Literal strings (not `DG.AGG` references) because the annotation
 *  generator only reads literal arrays — see the metadata traps in CLAUDE.md. */
export const AGGREGATION_TYPES = [
  'count', 'values', 'unique', 'nulls',
  'min', 'max', 'sum', 'avg', 'med', 'stdev', 'variance', 'q1', 'q2', 'q3', 'skew', 'kurt',
  'first', 'concat unique',
];

/** Aggregations that describe rows rather than a column's values — `count` is
 *  the only one that is meaningful with no column at all. */
const COLUMNLESS_AGGREGATIONS = new Set(['count']);

/** One entry of the `aggregations` JSON: which column, aggregated how, under
 *  what name. `column` may be blank only for {@link COLUMNLESS_AGGREGATIONS}. */
export interface AggregationSpec {
  column: string;
  type: string;
  name?: string;
}

// ---------- conditions ----------

/** Rows of `table` matching `condition`, as a bitset.
 *
 *  Evaluated by `core:AddNewColumn` with `aux.addColumn = false` — the flag its
 *  own `applyFormula` path uses to compute a formula column WITHOUT appending it
 *  (`ColumnFunctionsUtils.AddColumnAuxParam`, ddt/columns_functions.dart). So
 *  the user's table is never touched, and the expression grammar is exactly the
 *  one the formula editor autocompletes and validates. */
export async function conditionMask(table: DG.DataFrame, condition: string): Promise<DG.BitSet> {
  const expression = String(condition ?? '').trim();
  if (expression === '')
    throw new Error('Condition is empty — write a boolean expression, e.g. ${age} > 30');

  // `auto`, deliberately, NOT `bool`: asking for a bool COERCES — `${score} + 1`
  // comes back as a bool column of nulls and silently matches nothing. Letting
  // the engine infer the type is what makes the check below able to see that
  // the expression wasn't a condition at all.
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
  // `column.get(i)` is null for rows the formula could not evaluate; those are
  // "not matching", never "matching" — a null must not silently pass a filter.
  return DG.BitSet.create(table.rowCount, (i) => column.get(i) === true);
}

// ---------- rows ----------

/** Rows matching the condition, as a new table (core's `Subset`, row half). */
export async function filterRows(table: DG.DataFrame, condition: string): Promise<DG.DataFrame> {
  return named(table.clone(await conditionMask(table, condition)), `${table.name} filtered`);
}

/** Rows matching the condition removed — the complement of {@link filterRows}. */
export async function deleteRows(table: DG.DataFrame, condition: string): Promise<DG.DataFrame> {
  const mask = await conditionMask(table, condition);
  return named(table.clone(mask.invert()), `${table.name} trimmed`);
}

/** Rows matching the condition, keeping only the chosen columns — core's
 *  `Subset` in full. No columns chosen means every column. */
export async function extractRows(
  table: DG.DataFrame, condition: string, columns?: DG.Column[] | null,
): Promise<DG.DataFrame> {
  const names = columnNames(columns);
  const mask = await conditionMask(table, condition);
  return named(table.clone(mask, names.length > 0 ? names : null), `${table.name} extract`);
}

/** Selects the matching rows on the table and passes it on. Unlike the filters
 *  above this one IS a mutation — selection is table state, and copying the
 *  table to carry a selection would leave the selection on a copy nobody sees. */
export async function selectRows(
  table: DG.DataFrame, condition: string, clearSelection: boolean = true,
): Promise<DG.DataFrame> {
  const mask = await conditionMask(table, condition);
  if (clearSelection) table.selection.setAll(false, false);
  table.selection.or(mask);
  return table;
}

/** `count` random row indices out of `rowCount`, drawn without replacement from
 *  a seeded generator so a pipeline replays identically. Pure — unit-tested. */
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

/** A random sample of rows as a new table. */
export function filterRandomRows(table: DG.DataFrame, count: number, seed: number = 42): DG.DataFrame {
  const picked = new Set(randomIndices(table.rowCount, count, seed));
  return named(table.clone(DG.BitSet.create(table.rowCount, (i) => picked.has(i))),
    `${table.name} sample`);
}

/** Selects a random sample of rows on the table and passes it on. */
export function selectRandomRows(
  table: DG.DataFrame, count: number, seed: number = 42, clearSelection: boolean = true,
): DG.DataFrame {
  const picked = new Set(randomIndices(table.rowCount, count, seed));
  if (clearSelection) table.selection.setAll(false, false);
  for (const i of picked) table.selection.set(i, true, false);
  table.selection.fireChanged();
  return table;
}

// ---------- columns ----------

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

/** Sets a tag on the chosen columns and passes the table on. Tags are column
 *  metadata (units, `.formula`, a renderer hint), so there is nothing to copy —
 *  the mutation IS the result. */
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

/** Computes an expression into a new column of the table and returns it.
 *
 *  `core:AddNewColumn` does the same thing; the reason to have this is the
 *  parameter shape — the expression here carries Flow's inline formula editor
 *  (bound to the upstream table) instead of a bare text field, and the
 *  reactivity/error plumbing (`subscribeOnChanges`, `errorBehavior`) that a
 *  pipeline never wants is simply absent rather than hidden. */
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

// ---------- aggregate / unpivot ----------

/** Parse the `aggregations` JSON, tolerating a blank or half-typed value (this
 *  runs while the panel renders). Unknown aggregation types are dropped rather
 *  than thrown on — {@link aggregationProblems} is what reports them. Pure. */
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

/** What is wrong with an aggregation list, as user-facing labels (empty = fine).
 *  Shared by the node's readiness check and the function itself, so the hint and
 *  the runtime failure always agree. Pure and synchronous. */
export function aggregationProblems(specs: readonly AggregationSpec[]): string[] {
  if (specs.length === 0) return ['at least one aggregation'];
  const missing = specs
    .filter((a) => !COLUMNLESS_AGGREGATIONS.has(a.type) && a.column.trim() === '')
    .map((a) => `a column for "${a.type}"`);
  return Array.from(new Set(missing));
}

/** Group / pivot / aggregate — `core:Aggregate` with pickable inputs.
 *
 *  `aggregations` is a JSON list built by the aggregation editor; group-by and
 *  pivot are real column lists. Passing a pivot column is what turns this from a
 *  group-by into a pivot table, which is why the function is also findable as
 *  "pivot". */
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
    // Cast: the builder's TS signature is the `AGG` union, but the Dart side
    // takes the aggregation name as a string, which is how `concat unique`
    // (a `STR_AGG`, not an `AGG`) reaches it.
    builder.add(spec.type as DG.AggregationType, column, spec.name?.trim() || null);
  }
  for (const name of columnNames(pivotColumns))
    builder.pivot(name);
  return named(builder.aggregate(), `${table.name} aggregated`);
}

/** Wide → long: each `mergeColumns` column becomes a (category, value) row pair,
 *  with `copyColumns` repeated alongside. Delegates to `core:Unpivot` — the only
 *  thing wrong with that function on a canvas is that its column parameters are
 *  bare string lists, so this one takes real column slots. */
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

// ---------- helpers ----------

/** Names of a column slot's value, tolerating the null/undefined a blank
 *  `column_list` compiles to. */
function columnNames(columns: DG.Column[] | null | undefined): string[] {
  if (!columns) return [];
  return columns.filter((c) => c != null).map((c) => c.name);
}

/** Names a produced table, so output tabs and previews read as something other
 *  than `table`. */
function named(table: DG.DataFrame, name: string): DG.DataFrame {
  table.name = name;
  return table;
}

/** Deterministic PRNG (mulberry32) — `Math.random` would make every re-run of a
 *  sampling node produce a different table, which breaks both invalidation
 *  (a node re-runs and its downstream disagrees) and reproducibility. */
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
