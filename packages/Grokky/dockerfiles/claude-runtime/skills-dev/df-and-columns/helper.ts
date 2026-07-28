/**
 * Draft `grokky.*` helpers for the `datagrok-df-and-columns` skill.
 *
 * These will eventually move into `src/claude/grokky-api.ts`. Living here
 * keeps the skill review self-contained: SKILL.md cites these by name,
 * examples.ts type-checks against them, eval/prompts.json references them.
 */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

// TODO(grokky-api): `levenshtein` already exists at grokky-api.ts:100 but is not
// exported. Once we move these helpers in, drop this copy and reuse that one.
function levenshtein(a: string, b: string): number {
  if (a === b) return 0;
  if (!a.length) return b.length;
  if (!b.length) return a.length;
  const dp = Array.from({length: a.length + 1}, () => new Array<number>(b.length + 1).fill(0));
  for (let i = 0; i <= a.length; i++) dp[i][0] = i;
  for (let j = 0; j <= b.length; j++) dp[0][j] = j;
  for (let i = 1; i <= a.length; i++) {
    for (let j = 1; j <= b.length; j++) {
      const cost = a[i - 1] === b[j - 1] ? 0 : 1;
      dp[i][j] = Math.min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost);
    }
  }
  return dp[a.length][b.length];
}

// ---------------------------------------------------------------------------
// findColumns
// ---------------------------------------------------------------------------

export type FindColumnsQuery = {
  /** Substring or fuzzy target (case-insensitive). */
  name?: string;
  /** Single semType or array — any match counts. */
  semType?: string | string[];
  /** `DG.COLUMN_TYPE.*` or one of the broad bucket strings. */
  type?: DG.ColumnType | 'numerical' | 'categorical' | 'dateTime';
  /** Each entry: tag must have this exact value. `null`/`undefined` value matches any. */
  tag?: Record<string, string | null | undefined>;
  /** Enable Levenshtein fallback on names. Default true. */
  fuzzy?: boolean;
  /** Max candidates returned. Default 5. */
  limit?: number;
};

export type FindColumnsHit = {
  column: DG.Column;
  /** 0..1 confidence; 1 = exact name, ~0.5 = fuzzy. */
  score: number;
  /** Short human-readable explanation Claude can quote. */
  reason: string;
};

/**
 * Find columns matching a free-form intent. Tries (in order): exact name,
 * semType match, type bucket, tag predicate, fuzzy name match. Returns the
 * top `limit` ranked candidates. Column-discovery only — ignores filter and
 * selection.
 *
 * @example
 *   findColumns(df, {semType: DG.SEMTYPE.MOLECULE});
 *   findColumns(df, {name: 'mw', fuzzy: true});
 */
export function findColumns(df: DG.DataFrame, query: FindColumnsQuery): FindColumnsHit[] {
  const empty = !query.name && !query.semType && !query.type && !query.tag;
  if (empty) return [];
  const limit = query.limit ?? 5;
  const fuzzy = query.fuzzy !== false;
  const hits = new Map<DG.Column, FindColumnsHit>();
  const push = (col: DG.Column, score: number, reason: string) => {
    const prev = hits.get(col);
    if (!prev || prev.score < score) hits.set(col, {column: col, score, reason});
  };

  // 1. Type bucket — narrow the search space.
  let pool: DG.Column[];
  if (query.type === 'numerical') pool = df.columns.numerical as unknown as DG.Column[];
  else if (query.type === 'categorical') pool = df.columns.categorical as unknown as DG.Column[];
  else if (query.type === 'dateTime') pool = df.columns.dateTime as unknown as DG.Column[];
  else if (query.type) pool = df.columns.toList().filter((c) => c.type === query.type);
  else pool = df.columns.toList();

  // 2. Exact-name match (case-insensitive — same as df.col).
  if (query.name) {
    const wantedLc = query.name.toLowerCase();
    for (const c of pool) {
      if (c.name.toLowerCase() === wantedLc)
        push(c, 1.0, `exact name "${c.name}"`);
    }
  }

  // 3. semType match.
  if (query.semType) {
    const wanted = Array.isArray(query.semType) ? query.semType : [query.semType];
    for (const c of pool) {
      if (c.semType && wanted.includes(c.semType))
        push(c, 0.9, `semType match "${c.semType}"`);
    }
  }

  // 4. Tag predicate. `undefined`/`null` value matches any non-empty.
  if (query.tag) {
    const entries = Object.entries(query.tag);
    for (const c of pool) {
      let allMatch = true;
      const matched: string[] = [];
      for (const [k, v] of entries) {
        const has = c.tags.has(k);
        if (v === null || v === undefined) {
          if (!has) { allMatch = false; break; }
          matched.push(`${k}=*`);
        } else {
          if (c.getTag(k) !== v) { allMatch = false; break; }
          matched.push(`${k}=${v}`);
        }
      }
      if (allMatch && entries.length > 0)
        push(c, 0.8, `tag match ${matched.join(', ')}`);
    }
  }

  // 5. Substring match (case-insensitive).
  if (query.name) {
    const wantedLc = query.name.toLowerCase();
    for (const c of pool) {
      const nLc = c.name.toLowerCase();
      if (nLc === wantedLc) continue;     // already covered as exact
      if (nLc.includes(wantedLc))
        push(c, 0.7, `name contains "${query.name}"`);
    }
  }

  // 6. Fuzzy name match (Levenshtein).
  if (query.name && fuzzy) {
    const wantedLc = query.name.toLowerCase();
    const budget = Math.max(2, Math.ceil(wantedLc.length / 3));
    for (const c of pool) {
      const nLc = c.name.toLowerCase();
      if (nLc === wantedLc || nLc.includes(wantedLc)) continue;
      const d = levenshtein(nLc, wantedLc);
      if (d <= budget) {
        const score = 0.6 * (1 - d / (wantedLc.length + 1));
        push(c, Math.max(0.2, score), `fuzzy name distance ${d}`);
      }
    }
  }

  // 7. If only type/semType given and we have a pool, return everything from the type bucket
  //    as low-confidence candidates so Claude can still see the options.
  if (!query.name && !query.tag && hits.size === 0 && (query.type || query.semType)) {
    for (const c of pool)
      push(c, 0.5, query.type ? `type bucket "${query.type}"` : 'type bucket');
  }

  return Array.from(hits.values())
    .sort((a, b) => b.score - a.score)
    .slice(0, limit);
}

// ---------------------------------------------------------------------------
// describeColumn
// ---------------------------------------------------------------------------

export type ColumnDescription = {
  name: string;
  type: DG.ColumnType;
  semType: string | null;
  length: number;
  missing: number;
  unique: number;
  units?: string;
  format?: string;
  friendlyName?: string;
  description?: string;
  numerical?: {min: number; max: number; avg: number; stdev: number; med: number};
  categorical?: {topCategories: Array<{value: string; count: number}>};
};

/**
 * Concise JSON-friendly snapshot of a column. Uses `col.stats` which DG
 * caches per column. Does not cap on row count.
 *
 * @example
 *   describeColumn(df.getCol('activity'));
 */
export function describeColumn(col: DG.Column): ColumnDescription {
  const out: ColumnDescription = {
    name: col.name,
    type: col.type,
    semType: col.semType || null,
    length: col.length,
    missing: col.stats.missingValueCount,
    unique: col.stats.uniqueCount,
  };
  const units = col.meta.units; if (units) out.units = units;
  const format = col.meta.format; if (format) out.format = format;
  const friendlyName = col.meta.friendlyName; if (friendlyName) out.friendlyName = friendlyName;
  const description = col.meta.description; if (description) out.description = description;
  if (col.isNumerical) {
    const s = col.stats;
    out.numerical = {min: s.min, max: s.max, avg: s.avg, stdev: s.stdev, med: s.med};
  }
  if (col.type === DG.COLUMN_TYPE.STRING)
    out.categorical = {topCategories: topCategories(col, 5)};
  return out;
}

// ---------------------------------------------------------------------------
// topCategories
// ---------------------------------------------------------------------------

/**
 * Top-N category values with counts for a string column. Sorted by count desc.
 *
 * @example
 *   topCategories(df.getCol('country'), 10);
 */
export function topCategories(col: DG.Column, n: number = 5): Array<{value: string; count: number}> {
  if (col.type !== DG.COLUMN_TYPE.STRING) return [];
  const counts = new Map<string, number>();
  const len = col.length;
  for (let i = 0; i < len; i++) {
    const v = col.get(i) as string | null;
    if (v === null || v === undefined) continue;
    counts.set(v, (counts.get(v) ?? 0) + 1);
  }
  return Array.from(counts.entries())
    .sort((a, b) => b[1] - a[1])
    .slice(0, n)
    .map(([value, count]) => ({value, count}));
}

// ---------------------------------------------------------------------------
// setColumnMeta
// ---------------------------------------------------------------------------

export type ColorCodingSpec =
  | {kind: 'linear'; range?: string[]; min?: number; max?: number; belowMinColor?: string; aboveMaxColor?: string}
  | {kind: 'categorical'; map?: Record<string, string | number>; fallbackColor?: string | number}
  | {kind: 'conditional'; rules: Record<string, string | number>}
  | {kind: 'off'};

export type ColumnMetaSpec = {
  semType?: string | null;
  format?: string | null;
  units?: string | null;
  friendlyName?: string | null;
  description?: string | null;
  choices?: string[] | null;
  /** `undefined` = leave alone; `{kind:'off'}` = clear; anything else = apply. */
  colorCoding?: ColorCodingSpec;
  /** Raw tag escape hatch. `null` value deletes the tag. */
  tags?: Record<string, string | null>;
};

/**
 * Bulk-set column metadata. `null` clears a field; `undefined` leaves it alone.
 *
 * @example
 *   setColumnMeta(col, {semType: DG.SEMTYPE.Ki, units: 'nM', format: '0.00'});
 */
export function setColumnMeta(col: DG.Column, meta: ColumnMetaSpec): DG.Column {
  if (meta.semType !== undefined) col.semType = meta.semType as string;  // null also valid
  if (meta.format !== undefined) col.meta.format = meta.format;
  if (meta.units !== undefined) col.meta.units = meta.units;
  if (meta.friendlyName !== undefined) col.meta.friendlyName = meta.friendlyName;
  if (meta.description !== undefined) col.meta.description = meta.description;
  if (meta.choices !== undefined) col.meta.choices = meta.choices;

  if (meta.colorCoding !== undefined) {
    const cc = meta.colorCoding;
    if (cc.kind === 'linear') {
      if (!col.isNumerical)
        throw new Error(`setColumnMeta: linear color coding requires a numerical column; "${col.name}" is ${col.type}`);
      const opts: {min?: number; max?: number; belowMinColor?: string; aboveMaxColor?: string} = {};
      if (cc.min !== undefined) opts.min = cc.min;
      if (cc.max !== undefined) opts.max = cc.max;
      if (cc.belowMinColor !== undefined) opts.belowMinColor = cc.belowMinColor;
      if (cc.aboveMaxColor !== undefined) opts.aboveMaxColor = cc.aboveMaxColor;
      col.meta.colors.setLinear(cc.range ?? null, Object.keys(opts).length ? opts : null);
    } else if (cc.kind === 'categorical') {
      const opts = cc.fallbackColor !== undefined ? {fallbackColor: cc.fallbackColor} : null;
      col.meta.colors.setCategorical((cc.map ?? null) as any, opts as any);
    } else if (cc.kind === 'conditional') {
      col.meta.colors.setConditional(cc.rules as any);
    } else if (cc.kind === 'off') {
      col.meta.colors.setDisabled();
    }
  }

  if (meta.tags) {
    for (const [k, v] of Object.entries(meta.tags)) {
      if (v === null) col.tags.delete(k);
      else col.setTag(k, v);
    }
  }
  return col;
}

// ---------------------------------------------------------------------------
// cloneDf
// ---------------------------------------------------------------------------

export type CloneDfOpts = {
  rows?: DG.BitSet | 'filtered' | 'selected';
  cols?: string[];
  withSelection?: boolean;
  withTags?: boolean;
};

/**
 * Pure rewrite of `df.clone(...)` with named args. `rows: 'filtered'` uses
 * `df.filter`; `'selected'` uses `df.selection`.
 *
 * @example
 *   cloneDf(df, {cols: ['smiles', 'activity'], rows: 'filtered'});
 */
export function cloneDf(df: DG.DataFrame, opts: CloneDfOpts = {}): DG.DataFrame {
  let mask: DG.BitSet | null = null;
  if (opts.rows === 'filtered') mask = df.filter;
  else if (opts.rows === 'selected') mask = df.selection;
  else if (opts.rows instanceof DG.BitSet) mask = opts.rows;
  return df.clone(mask, opts.cols ?? null, opts.withSelection ?? false, opts.withTags ?? true);
}

// ---------------------------------------------------------------------------
// addColumn
// ---------------------------------------------------------------------------

export type AddColumnSpec = {
  name: string;
  type?: DG.ColumnType;
  formula?: string;
  values?: any[];
  init?: (i: number) => any;
  virtual?: (i: number) => any;
  insertAt?: number;
  ensureUniqueName?: boolean;
  meta?: ColumnMetaSpec;
};

/**
 * Add a column choosing the right backing call. Always async for uniformity.
 * At most one of `formula`/`values`/`init`/`virtual` may be set.
 *
 * @example
 *   await addColumn(df, {name: 'Ki', type: DG.COLUMN_TYPE.FLOAT,
 *     meta: {semType: DG.SEMTYPE.Ki, units: 'nM'}});
 */
export async function addColumn(df: DG.DataFrame, spec: AddColumnSpec): Promise<DG.Column> {
  const sources = [spec.formula, spec.values, spec.init, spec.virtual].filter((s) => s !== undefined);
  if (sources.length > 1)
    throw new Error('addColumn: only one of formula / values / init / virtual may be set');

  const name = (spec.ensureUniqueName ?? true) ? df.columns.getUnusedName(spec.name) : spec.name;
  let col: DG.Column;

  if (spec.formula !== undefined) {
    col = await df.columns.addNewCalculated(name, spec.formula, spec.type ?? 'auto');
  } else if (spec.values !== undefined) {
    const t = spec.type ?? inferTypeFromValues(spec.values);
    const built = DG.Column.fromList(t, name, spec.values);
    col = df.columns.add(built);
  } else if (spec.virtual !== undefined) {
    col = df.columns.addNewVirtual(name, spec.virtual, spec.type ?? DG.COLUMN_TYPE.OBJECT);
  } else if (spec.init !== undefined) {
    col = df.columns.addNew(name, spec.type ?? DG.COLUMN_TYPE.FLOAT);
    col.init(spec.init);
  } else {
    col = df.columns.addNew(name, spec.type ?? DG.COLUMN_TYPE.STRING);
  }

  if (spec.insertAt !== undefined && spec.insertAt >= 0) {
    // addNew/add appended at end — re-insert at the wanted position.
    df.columns.remove(col);
    df.columns.insert(col, spec.insertAt);
  }

  if (spec.meta) setColumnMeta(col, spec.meta);
  return Promise.resolve(col);
}

function inferTypeFromValues(values: any[]): DG.ColumnType {
  for (const v of values) {
    if (v === null || v === undefined) continue;
    if (typeof v === 'number') return Number.isInteger(v) ? DG.COLUMN_TYPE.INT : DG.COLUMN_TYPE.FLOAT;
    if (typeof v === 'boolean') return DG.COLUMN_TYPE.BOOL;
    if (v instanceof Date) return DG.COLUMN_TYPE.DATE_TIME;
    return DG.COLUMN_TYPE.STRING;
  }
  return DG.COLUMN_TYPE.STRING;
}

// ---------------------------------------------------------------------------
// removeColumns
// ---------------------------------------------------------------------------

/**
 * Bulk remove columns by name. `onMissing` controls what happens when a name
 * is not present: `'skip'` (default) silently ignores; `'warn'` logs;
 * `'throw'` raises.
 *
 * @returns The list of names that were actually removed.
 *
 * @example
 *   removeColumns(df, ['_tmp1', '_tmp2']);
 */
export function removeColumns(
  df: DG.DataFrame,
  names: string[],
  onMissing: 'skip' | 'throw' | 'warn' = 'skip',
): string[] {
  const removed: string[] = [];
  for (const name of names) {
    const col = df.col(name);
    if (!col) {
      if (onMissing === 'throw') throw new Error(`removeColumns: column "${name}" not found`);
      if (onMissing === 'warn') console.warn(`removeColumns: column "${name}" not found, skipping`);
      continue;
    }
    df.columns.remove(col);
    removed.push(col.name);
  }
  return removed;
}

// ---------------------------------------------------------------------------
// renameColumn
// ---------------------------------------------------------------------------

/**
 * Rename a column. With `ensureUnique`, runs through `getUnusedName` to avoid
 * collisions and returns the final name actually used. Warns about viewer
 * bindings — does not update them.
 *
 * @example
 *   renameColumn(df, 'pIC50', 'pIC50_old', {ensureUnique: true});
 */
export function renameColumn(
  df: DG.DataFrame,
  from: string,
  to: string,
  opts: {ensureUnique?: boolean} = {},
): string {
  const col = df.getCol(from);  // throws if missing — fail-fast intent
  const finalName = opts.ensureUnique ? df.columns.getUnusedName(to) : to;
  if (finalName === col.name) return finalName;
  // Detect viewer bindings (best-effort: scan TableView viewers for column refs by old name).
  try {
    const tv = grok.shell.tv;
    if (tv && tv.dataFrame === df) {
      const bound: string[] = [];
      for (const v of tv.viewers) {
        const opts = (v as any).getOptions?.()?.look ?? {};
        for (const [k, val] of Object.entries(opts)) {
          if (typeof val === 'string' && val === from) bound.push(`${v.type}.${k}`);
        }
      }
      if (bound.length > 0)
        console.warn(`renameColumn: "${from}" is referenced by viewer property bindings (${bound.join(', ')}). These will not be auto-updated.`);
    }
  } catch { /* viewer inspection is best-effort */ }
  col.name = finalName;
  return finalName;
}
