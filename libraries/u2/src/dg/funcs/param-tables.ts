/* The per-table-param machinery behind FuncCallForm's W3 routes (dataframe / column /
   column_list): the dependent retargeting on the param stream, the Dart default-column pick
   (fpe:162-187), the `ColumnFilter.fromProp`-parity predicate and the feature-detected metadata
   reader. Internal to `src/dg/funcs` — func-form.ts is the only importer; not exported through
   `src/dg/index.ts`. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {batch} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {span} from '../../core/elements.js';
import {Input} from '../../core/input-base.js';
import {ColumnInput} from '../inputs/column-combo.js';
import {ColumnsInput} from '../inputs/columns.js';
import {ObjectForm} from '../forms/object-form.js';
import type {FuncCallLike, FuncCallParamLike, FuncField} from './func-form.js';

const api = globalThis as {grok_Property_Get?: (dart: unknown, name: string) => unknown};

/** Platform metadata js-api has no wrapper getter for, through the generic registration the
 * client installs (`grok_api.dart:898`) — the `grok_Property_Get_PropertySubType` precedent.
 * Null when the global or the handle is absent — and on an undeclared key, where the Dart
 * `getPropertyValue` THROWS (throwOnError defaults true) rather than answering null. */
export function propMeta(prop: {dart?: unknown}, name: string): unknown {
  const dart = prop?.dart;
  if (dart == null || typeof api.grok_Property_Get !== 'function')
    return null;
  try {
    return api.grok_Property_Get(dart, name) ?? null;
  } catch {
    return null;
  }
}

/** The annotation linker's field first — the single source of truth: the implicit link (`column
 * age` after `dataframe df`) never reaches `options['table']` — then the explicit option, which
 * also serves signature-registered funcs decorated at runtime. */
export function parentTableName(prop: FuncCallParamLike['property']): string | null {
  const linked = propMeta(prop, 'parentTableParamName');
  if (typeof linked === 'string' && linked !== '')
    return linked;
  const option = prop.options?.['table'];
  return typeof option === 'string' && option !== '' ? option : null;
}

/** A DataFrame/Column readback, or null for anything else: a Resolve\* FuncCall — a raw string
 * written by external code, pending resolution — must not leak a bogus name (P-W3-2). */
export function asNamed(v: any): {name: string} | null {
  return v != null && typeof v.name === 'string' && !('func' in v) ? v : null;
}

export function asTable(v: any): DG.DataFrame | null {
  return asNamed(v) != null && v.columns != null ? v as DG.DataFrame : null;
}

/** A param value as the frame it denotes: a name (a custom table editor may write one) resolves
 * through the shell, anything else must already be a frame. */
export function resolveTable(v: any): DG.DataFrame | null {
  return typeof v === 'string' ? tableByName(v) : asTable(v);
}

function openTables(): DG.DataFrame[] {
  return grok.shell.tables ?? [];
}

/** Name → the open table (shell names are unique — `addTable` dedupes); one interop call,
 * where scanning `grok.shell.tables` wraps every open frame. */
export function tableByName(name: string | null): DG.DataFrame | null {
  return name == null ? null : (grok.shell.tableByName(name) as DG.DataFrame | null) ?? null;
}

/** The auto-fill source (fpe:62-64): the current table, else the first open one. */
export function preferredTable(): DG.DataFrame | null {
  return (grok.shell.t as DG.DataFrame | null) ?? openTables()[0] ?? null;
}

/** A `column_list` readback as names: a ColumnList (P-W3-9), an array of columns, or any
 * degenerate state. */
export function columnNamesOf(v: any): string[] {
  if (v == null)
    return [];
  if (typeof v.names === 'function')
    return v.names();
  return Array.isArray(v) ?
    v.map((c) => asNamed(c)?.name).filter((n): n is string => n != null) : [];
}

const HAS_MISSING = 'has missing values';
const NO_MISSING = 'no missing values';

/** `ColumnFilter.fromProp` parity (`column_filter.dart:32-82`): the type/magic-string filter AND
 * the semType match AND the maxCategories cap (read through the metadata door). */
export function columnPredicate(prop: FuncCallParamLike['property']): (c: DG.Column) => boolean {
  const semType = prop.semType;
  const rawMax = propMeta(prop, 'maxCategories');
  const max = typeof rawMax === 'number' && isFinite(rawMax) ? rawMax : null;
  const byType = typeFilter(prop.columnTypeFilter ?? null);
  return (c) => byType(c) &&
    (semType == null || semType === '' || c.semType === semType) &&
    (c.isCategorical !== true || max === null || (c.categories?.length ?? 0) <= max);
}

function typeFilter(filter: string | null): (c: DG.Column) => boolean {
  if (filter == null || filter === '')
    return () => true;
  if (filter === HAS_MISSING || filter === NO_MISSING) {
    const expected = filter === HAS_MISSING;
    return (c) => ((c.stats?.missingValueCount ?? 0) > 0) === expected;
  }
  switch (filter) {
    case 'numerical':
      return (c) => c.isNumerical === true;
    case 'numerical_no_datetime':
      return (c) => c.isNumerical === true && c.type !== 'datetime';
    case 'categorical':
      return (c) => c.isCategorical === true;
    case 'datetime':
      return (c) => c.type === 'datetime';
    case 'categorical_or_datetime':
      return (c) => c.isCategorical === true || c.type === 'datetime';
    default:
      return (c) => c.type === filter;
  }
}

const KIND_WORDS: Record<string, string> = {
  numerical: 'numerical', numerical_no_datetime: 'numerical', categorical: 'categorical',
  datetime: 'datetime', categorical_or_datetime: 'categorical or datetime',
};

/** Why a column field has nothing to offer — the filter kind (or semType) and the table it
 * scanned, in the user's words. */
export function noMatchMessage(prop: FuncCallParamLike['property'], tableName: string): string {
  const semType = prop.semType;
  if (semType != null && semType !== '')
    return `No ${semType} columns in '${tableName}'`;
  const filter: string | null = prop.columnTypeFilter ?? null;
  if (filter == null || filter === '')
    return `No columns in '${tableName}'`;
  if (filter === HAS_MISSING || filter === NO_MISSING)
    return `No columns ${filter === HAS_MISSING ? 'with' : 'without'} missing values in '${tableName}'`;
  return `No ${KIND_WORDS[filter] ?? filter} columns in '${tableName}'`;
}

const MOLECULE_NAMES = new Set(['structure', 'smiles', 'canonical_smiles', 'molecule']);

/** Structural, like column-combo's platform door: a vendored datagrok-api may predate the 1.28
 * `StringUtils` distance statics. */
const strings = DG.StringUtils as unknown as {
  levenshteinDistance(a: string, b: string): number,
  jaroWinklerDistance(a: string, b: string): number,
};

/** The Dart default pick (fpe:162-179, exact): among the passing columns — none → null; a
 * semType-carrying param takes the molecule-named or first one; otherwise the nearest lowercased
 * name by normalized distance — Levenshtein for 1-character param names, Jaro-Winkler
 * otherwise. */
export function pickColumn(table: DG.DataFrame, paramName: string,
  prop: FuncCallParamLike['property'], filter: (c: DG.Column) => boolean): DG.Column | null {
  const candidates = table.columns.toList().filter(filter);
  if (candidates.length === 0)
    return null;
  if (prop.semType != null && prop.semType !== '')
    return candidates.find((c) => MOLECULE_NAMES.has(c.name.toLowerCase())) ?? candidates[0];
  const name = paramName.toLowerCase();
  const distance = name.length === 1 ? strings.levenshteinDistance : strings.jaroWinklerDistance;
  let best = candidates[0];
  let bestDistance = distance(name, best.name.toLowerCase());
  for (let i = 1; i < candidates.length; i++) {
    const d = distance(name, candidates[i].name.toLowerCase());
    if (d < bestDistance) {
      best = candidates[i];
      bestDistance = d;
    }
  }
  return best;
}

export interface ColumnDependent {
  kind: 'column' | 'columns';
  field: FuncField;
  filter: (c: DG.Column) => boolean;
}

/** The auto affordance: a value written by an auto-pick or the table auto-fill is a guess, and
 * wears a dimmed `auto` badge until the user edits the field themselves — a DOM `change` only a
 * real interaction fires; a binding or system write never does, so the guess mark survives
 * retargets. The first user edit clears it for good. */
export function markAuto(field: FuncField): void {
  if (field.userTouched === true)
    return;
  if (field.autoBadge === undefined) {
    const badge = span('auto', 'u2-param-auto');
    badge.title = 'Chosen automatically';
    field.autoBadge = badge;
    const touched = () => {
      field.userTouched = true;
      badge.remove();
      field.input.root.removeEventListener('change', touched);
    };
    field.input.root.addEventListener('change', touched);
  }
  if (field.autoBadge.parentElement === null)
    field.input.box.append(field.autoBadge);
}

export function unmarkAuto(field: FuncField): void {
  field.autoBadge?.remove();
}

/** The dependent-rewiring owner: one instance per rendered table param, constructed under the
 * form's scope after all fields exist, disposed and rebuilt on a `source` rebind. One
 * subscription on the PARAM stream (FE-7 #2/#3 — divergence #12: both a user pick and an
 * external `setParamValue` retarget; Dart wires the input and stacks a listener per dependent),
 * so same-frame writes cannot echo — they are suppressed at the source (P-W3-2). */
export class TableBinding {
  private readonly _call: FuncCallLike;
  private readonly _tableParam: FuncCallParamLike;
  private readonly _tableOf: (name: string | null) => DG.DataFrame | null;
  private readonly _deps: ColumnDependent[];
  private readonly _notify: (dep: ColumnDependent, oldTableName: string) => void;
  private readonly _disown: () => void;
  private _tableSub: {unsubscribe(): void} | undefined;
  private _disposed = false;

  constructor(call: FuncCallLike, tableParam: FuncCallParamLike,
    tableOf: (name: string | null) => DG.DataFrame | null, dependents: ColumnDependent[],
    notify: (dep: ColumnDependent, oldTableName: string) => void) {
    this._call = call;
    this._tableParam = tableParam;
    this._tableOf = tableOf;
    this._deps = dependents;
    this._notify = notify;
    const owner = Scope.ambient!;
    const cleanup = () => this._dispose();
    owner.own(cleanup);
    this._disown = () => owner.disown(cleanup);
    const sub = tableParam.onChanged.subscribe(() => this._retargetAll());
    this._tableSub = sub;
  }

  /** Initial retarget to the param's current table + the column auto-pick (fpe:185-187): a
   * null-valued `column` dependent gets the pick WRITTEN into the call (`??=` in Dart is a real
   * write, P-W3-3; unguarded by `skipDefaultInit` — Dart's guard covers only
   * `options['default']`, fpe:558); column lists are never auto-picked. A pre-set `columns`
   * value survives: `changeTable` runs only when the bound table differs, and the param's names
   * are loaded into the input otherwise (P-W3-9). */
  start(): void {
    const table = this._table();
    for (const dep of this._deps) {
      if (dep.kind === 'columns')
        this._startColumns(dep, table);
      else
        this._startColumn(dep, table);
      // the required verdict reads the table, which can move under an unchanged (null) value
      dep.field.input.revalidate();
    }
  }

  dispose(): void {
    this._disown();
    this._dispose();
  }

  private _table(): DG.DataFrame | null {
    const v = this._tableParam.value;
    return typeof v === 'string' ? this._tableOf(v) : asTable(v);
  }

  /** js-api wrappers are re-minted per access, so wrapper identity is meaningless; the dart
   * handle (or, failing that, the shell-unique name) is the frame's stable identity. */
  private static _sameTable(a: DG.DataFrame | null, b: DG.DataFrame | null): boolean {
    if (a === b)
      return true;
    if (a == null || b == null)
      return false;
    const da = a.dart;
    const db = b.dart;
    return da != null && db != null ? da === db : a.name === b.name;
  }

  /** The control follows its frame itself (ColumnInput's internal subscriptions) — the binding
   * only retargets it and re-defaults, mirroring {@link _startColumns}. */
  private _startColumn(dep: ColumnDependent, table: DG.DataFrame | null): void {
    const input = dep.field.input as ColumnInput;
    batch(() => {
      if (!TableBinding._sameTable(input.table, table))
        input.changeTable(table, dep.filter);
      if (dep.field.param.value != null || table === null) {
        const name = asNamed(dep.field.param.value)?.name ?? null;
        if (name !== input.value.peek())
          Input.system(() => input.value.value = name);
        unmarkAuto(dep.field);
        return;
      }
      const picked = pickColumn(table, dep.field.param.name, dep.field.param.property, dep.filter);
      if (picked === null) {
        unmarkAuto(dep.field);
        return;
      }
      // the call first: the display write then reads back its own name and does not echo
      this._call.setParamValue(dep.field.param.name, picked);
      Input.system(() => input.value.value = picked.name);
      markAuto(dep.field);
    });
  }

  private _startColumns(dep: ColumnDependent, table: DG.DataFrame | null): void {
    const input = dep.field.input as ColumnsInput;
    batch(() => {
      if (!TableBinding._sameTable(input.table, table))
        input.changeTable(table, dep.filter);
      const names = columnNamesOf(dep.field.param.value);
      if (!ObjectForm.same(names, input.value.peek()))
        Input.system(() => input.value.value = names);
    });
  }

  /** Table change → every dependent retargets and re-defaults, writing into the call —
   * overwrite incl. null on no-match is Dart's own behavior (P-W3-4, P-W3-10). */
  private _retargetAll(): void {
    const table = this._table();
    for (const dep of this._deps) {
      if (dep.kind === 'columns')
        this._retargetColumns(dep, table);
      else
        this._retargetColumn(dep, table);
      dep.field.input.revalidate();
    }
  }

  /** The `[]` write is Dart's own behavior (`columns_input.dart:118-129`) — but a non-empty
   * selection must not vanish in silence: the field says whose columns it lost. */
  private _retargetColumns(dep: ColumnDependent, table: DG.DataFrame | null): void {
    const input = dep.field.input as ColumnsInput;
    const cleared = input.value.peek();
    const from = input.table;
    input.changeTable(table, dep.filter);
    if (cleared.length > 0 && from != null)
      this._notify(dep, from.name);
  }

  private _retargetColumn(dep: ColumnDependent, table: DG.DataFrame | null): void {
    const input = dep.field.input as ColumnInput;
    const picked = table === null ? null :
      pickColumn(table, dep.field.param.name, dep.field.param.property, dep.filter);
    const before = input.value.peek();
    // one batch, so changeTable's internal clear and the re-default flush as ONE field-effect
    // write — the re-default, not the clear's null, is what reaches the call
    batch(() => {
      input.changeTable(table, dep.filter);
      input.value.value = picked?.name ?? null;
    });
    // an unchanged name is still a new column: the effect saw no edit, and the call must not
    // keep the old table's object
    if (picked !== null && input.value.peek() === before)
      this._call.setParamValue(dep.field.param.name, picked);
    if (picked !== null)
      markAuto(dep.field);
    else
      unmarkAuto(dep.field);
  }

  private _dispose(): void {
    if (this._disposed)
      return;
    this._disposed = true;
    this._tableSub?.unsubscribe();
    this._tableSub = undefined;
  }
}
