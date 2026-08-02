import type {Dayjs} from 'dayjs';

/** Column reference: a declared/system column of the table, or a dotted FK path
 * ('project_id.name', up to 3 forward hops). */
export type DomainColumnRef<TColumn extends string = string> = TColumn | `${string}.${string}`;

/** Operators of the canonical condition tree (server: filter_compiler.dart:149-153 + fuzzy). */
export type DomainConditionOperator =
  '=' | '!=' | '>' | '>=' | '<' | '<=' | 'like' | 'not like' | '~*' | '!~*' | 'is' | 'is not' | 'fuzzy';

export type DomainFilterValue =
  string | number | boolean | null | Dayjs | (string | number | boolean)[];

/** One condition node. `value` may be a list (= ANY / != ALL), null (IS NULL / IS NOT NULL
 * under '='/'is'), or '@current' for the current user on user columns. Values are ALWAYS
 * bound server-side — never interpolate values into filter strings. */
export interface DomainCondition<TColumn extends string = string> {
  property: DomainColumnRef<TColumn>;
  operator: DomainConditionOperator;
  value?: DomainFilterValue;
}

/** Inside a tree, a bare string is a connector; an omitted connector means 'and'
 * (sticky — the last seen connector repeats). */
export type DomainConditionNode<TColumn extends string = string> =
  DomainCondition<TColumn> | 'and' | 'or' | DomainConditionTree<TColumn>;
export type DomainConditionTree<TColumn extends string = string> =
  DomainConditionNode<TColumn>[];

/** Every filter slot accepts the smart-filter string (human shorthand), a single condition,
 * or the canonical condition tree — uniformly across query, queryDf, aggregate, facets,
 * count, exists, and deleteWhere (server: DomainFilterCompiler.toTree). */
export type DomainFilter<TColumn extends string = string> =
  string | DomainCondition<TColumn> | DomainConditionTree<TColumn>;

export const DOMAIN_SYSTEM_COLUMNS = ['id', 'version', 'created_on', 'updated_on', 'author_id'] as const;

export interface DomainColumnError { column: string; code: string; message: string; }
export interface DomainRowErrors { index: number; id?: string | null; errors: DomainColumnError[]; }

export interface DomainInsertResult {
  id: string;
  created: boolean;
  version?: number;
  /** Set when the row was NOT created: 'duplicate' (business-key match) or 'idempotent-replay'. */
  status?: 'duplicate' | 'idempotent-replay';
  existingId?: string;
}
export interface DomainUpdateResult { id: string; version: number; }
export interface DomainDeleteResult { id: string; deleted: true; }
export type DomainOpResult = DomainInsertResult | DomainUpdateResult | DomainDeleteResult;

/** Wire shape of _audit rows (repository.dart rowAudit/tableAudit/schemaAudit selects). */
export interface DomainAuditEntry {
  /** Audit sequence number (domain_audit.id bigint — a wire number, not a uuid). */
  id: number;
  /** PG transaction id grouping multi-op writes (bigint). */
  tx_id: number | null;
  row_id?: string | null;    // present on table- and schema-wide audits
  table_id?: string | null;  // present on schema audits ('ddl' rows carry null)
  op: 'insert' | 'update' | 'delete' | 'promote' | 'ddl';
  actor_id: string | null;
  session_id: string | null;
  source: string;
  ts: string;                // ISO-8601
  /** Column-value map before the write; null on 'insert' ops (typed `any` so existing
   * unguarded consumers compile — guard before dereferencing). */
  before: any;
  /** Column-value map after the write; null on 'delete' ops (typed `any` — see `before`). */
  after: any;
}

/** Report of domain table `deleteWhere`: rows actually soft-deleted in this call,
 * and whether more matching deletable rows remain (loop while `hasMore`). */
export interface DomainDeleteReport { deleted: number; hasMore: boolean; }

export type DomainPermission = 'View' | 'Edit' | 'Delete' | 'Share';

/** One direct permission row on a domain registry entity (see `DomainTableClient.grants`). */
export interface DomainGrant {
  group: {id: string; friendlyName: string; personal: boolean};
  permission: DomainPermission;
}

export interface DomainBatchRowResult {
  index: number;
  id: string | null;
  status: 'inserted' | 'updated' | 'duplicate' | 'error';
  existingId?: string;
  errors?: DomainColumnError[];
}

/** Base of all structured /domains failures; [body] is the server's error envelope. */
export class DomainError extends Error {
  constructor(message: string, public readonly status: number,
              public readonly body: {[key: string]: any}) {
    super(message);
    this.name = new.target.name;
  }
  /** Server discriminant: 'validation' | 'version-conflict' | 'restrict' | 'filter' |
   * 'forbidden' | 'not-found' | 'manifest-validation' | 'invalid-mode' | 'id-collision' |
   * 'destructive-confirmation-required' | schema-mgmt codes | '' (transport). */
  get code(): string { return `${this.body['error'] ?? ''}`; }
  /** Index of the failing op in a transaction ops list; undefined otherwise. */
  get opIndex(): number | undefined { return this.body['opIndex']; }
}

export class DomainValidationError extends DomainError {          // code 'validation'
  get rows(): DomainRowErrors[] { return this.body['rows'] ?? []; }
  /** True when the failure is a business-key/idempotency uniqueness conflict
   * (per-row error code 'unique' — how the server reports duplicates under errorOnDuplicate). */
  get isDuplicate(): boolean {
    return this.rows.some((r) => (r.errors ?? []).some((e) => e.code === 'unique'));
  }
}
export class DomainVersionConflictError extends DomainError {     // code 'version-conflict'
  get id(): string { return this.body['id']; }
  get currentVersion(): number { return this.body['currentVersion']; }
  get expectedVersion(): number { return this.body['expectedVersion']; }
}
export class DomainRestrictError extends DomainError {}           // code 'restrict'
export class DomainFilterError extends DomainError {}             // code 'filter'
export class DomainForbiddenError extends DomainError {}          // code 'forbidden'
export class DomainNotFoundError extends DomainError {}           // code 'not-found'
export class DomainManifestValidationError extends DomainError {  // code 'manifest-validation'
  get errors(): {[key: string]: any}[] { return this.body['errors'] ?? []; }
}

/** Maps the interop's '#domainError' plain object to a typed class; passes anything else through. */
export function toDomainError(e: any): any {
  if (e == null || e['#domainError'] !== true)
    return e;
  const body = e.body ?? {};
  const ctor = ({
    'validation': DomainValidationError,
    'version-conflict': DomainVersionConflictError,
    'restrict': DomainRestrictError,
    'filter': DomainFilterError,
    'forbidden': DomainForbiddenError,
    'not-found': DomainNotFoundError,
    'manifest-validation': DomainManifestValidationError,
  } as {[code: string]: typeof DomainError})[`${body['error']}`] ?? DomainError;
  return new ctor(`${e.message ?? body['message'] ?? 'Domain request failed'}`, e.status ?? 0, body);
}

/** Awaits [p], rethrowing interop-marked failures as typed DomainErrors. */
export async function domainCall<T>(p: Promise<T>): Promise<T> {
  try {
    return await p;
  } catch (e: any) {
    throw toDomainError(e);
  }
}

/** Query options for domain table `query`/`queryDf`. */
export interface DomainQuerySpec<TColumn extends string = string, TExpandKey extends string = string> {
  /** Smart-filter string (same grammar as entity search, e.g. `barcode starts "P-1"`),
   * a single condition, or the canonical condition tree — values are bound server-side. */
  filter?: DomainFilter<TColumn>;
  /** Comma-separated column list; `!` prefix for descending, e.g. `'name,!created_on'`. */
  sort?: string;
  /** Columns to return; omit for all viewable columns. */
  columns?: TColumn[];
  /** Expansions (same-schema, depth 1): `'<fk_column>'` returns the master row's declared
   * columns prefixed `'<fk_column>.<name>'`; `'details:<table>[.<fk_column>]'` returns capped
   * child-row arrays under the child-table name (JSON queries only — not supported by
   * `queryDf`). */
  expand?: TExpandKey[];
  limit?: number;
  offset?: number;
}

/** One measure of {@link DomainAggregateSpec}: `fn` over `column` (`count` needs no column);
 * the output name is `as` (defaults to `fn` or `<fn>_<column>`). */
export interface DomainAggregateMeasure<TColumn extends string = string, TAlias extends string = string> {
  fn: 'count' | 'sum' | 'avg' | 'min' | 'max';
  column?: DomainColumnRef<TColumn>;
  as?: TAlias;
}

/** Spec for domain table `aggregate`. */
export interface DomainAggregateSpec<TColumn extends string = string,
    TGroup extends string = string, TAlias extends string = string> {
  /** Column names to group by; omit for a single-row grand total. */
  groupBy?: (DomainColumnRef<TColumn> & TGroup)[];
  measures: DomainAggregateMeasure<TColumn, TAlias>[];
  /** Filter applied before aggregation (smart string, condition, or condition tree). */
  filter?: DomainFilter<TColumn>;
  /** Comma-separated output names (group columns or measure aliases); `!` prefix for descending. */
  sort?: string;
  limit?: number;
  /** Master-FK expands ('<fk_column>'); groupBy/measures may then use '<fk>.<col>' (§13.2). */
  expand?: string[];
}
export type DomainAggregateRow<TKeys extends string = string> =
  {[K in TKeys]: number | string | boolean | null};

export type DomainFacetKind = 'categories' | 'histogram' | 'minMax' | 'count' | 'plan';

/** One facet request of {@link DomainFacetsSpec}; `id` keys its result in the response. */
export interface DomainFacetSpec<TColumn extends string = string,
    TId extends string = string, TKind extends DomainFacetKind = DomainFacetKind> {
  /** Response key for this facet's result. */
  id: TId;
  kind: TKind;
  /** Column name; `'categories'` also accepts a dotted FK path (e.g. `'category_id.name'`, up to
   * 3 hops — the counts respect the referenced table's row predicate). Not used by `'count'`/`'plan'`. */
  column?: DomainColumnRef<TColumn>;
  /** `'plan'` only: columns to profile (capped distinct count plus numeric/datetime min/max). */
  columns?: TColumn[];
  /** `'categories'` only: category cap (default 100, clamped to 1..1000); the result carries
   * `hasMore` when more remain — narrow with `search` instead of raising the cap. */
  limit?: number;
  /** `'categories'` only: server-side substring narrowing (compiled as a bound ILIKE). */
  search?: string;
  /** `'histogram'` only: bucket count (default 20, clamped to 1..200). */
  bins?: number;
  /** `'histogram'` only: pinned lower bound (a number, or an ISO-8601 string for datetime columns)
   * so zooming does not re-derive the axis; omit to bucket over the data bounds. */
  min?: number | string;
  /** `'histogram'` only: pinned upper bound (see `min`). */
  max?: number | string;
}

/** Spec for domain table `facets`: one request computes every facet in a single round
 * trip. Category counts, histogram buckets, and `'count'` are computed under `filter` with the
 * conditions on that facet's own column stripped, so a filter control shows counts under all
 * OTHER filters (classic faceted search). Exception: `'minMax'`, `'plan'`, and histogram BOUNDS
 * are computed under the row predicate only — the stable-axis rule ignores `filter` so a
 * narrowing filter never re-derives the axis. At most 32 facets per request. */
export interface DomainFacetsSpec<TColumn extends string = string,
    TId extends string = string, TKind extends DomainFacetKind = DomainFacetKind> {
  /** Smart string, single condition, or condition tree; omit for unfiltered counts. */
  filter?: DomainFilter<TColumn>;
  facets: DomainFacetSpec<TColumn, TId, TKind>[];
}

/** One category bucket of a `'categories'` facet: `total` ignores the filter, `filtered` respects
 * it (minus the facet's own column); ref columns group by id and carry the referenced row's
 * display name in `display`. */
export interface DomainFacetCategory { value: any; display?: string; total: number; filtered: number; }
export interface DomainFacetCategoriesResult { categories: DomainFacetCategory[]; hasMore?: boolean; }
export interface DomainFacetHistogramResult {
  min: number | string | null; max: number | string | null;
  buckets: number[]; totalBuckets: number[]; nulls: number;
}
export interface DomainFacetMinMaxResult { min: number | string | null; max: number | string | null; }
export interface DomainFacetCountResult { count: number; }
export interface DomainFacetPlanResult {
  columns: {name: string; distinct: number; min?: number | string; max?: number | string}[];
}
export type DomainFacetResultOf<K extends DomainFacetKind> =
  K extends 'categories' ? DomainFacetCategoriesResult :
  K extends 'histogram' ? DomainFacetHistogramResult :
  K extends 'minMax' ? DomainFacetMinMaxResult :
  K extends 'count' ? DomainFacetCountResult : DomainFacetPlanResult;

/** One operation of `DomainsDataSource.transaction`. */
export interface DomainTransactionOp {
  op: 'insert' | 'update' | 'delete';
  /** Table name within the transaction's schema. */
  table: string;
  /** Names this op's new row id; later ops may reference it in values as `'$<ref>'`
   * (escape a literal leading `$` in a value by doubling it: `'$$100'` stores `'$100'`). */
  ref?: string;
  values?: object;
  id?: string;
  /** Optimistic-concurrency guard for update ops. */
  expectedVersion?: number;
}

/** Options for domain table `batch`. */
export interface DomainBatchOptions {
  /** `'insert'` (default) or `'upsert'` (merge by the table's business key). */
  mode?: 'insert' | 'upsert';
  /** Abort the whole batch on any row error (default true); false applies good rows
   * and reports bad ones per row. */
  allOrNothing?: boolean;
  /** Report business-key duplicates as errors instead of skipping them. */
  errorOnDuplicate?: boolean;
  /** Payload format for `Uint8Array` data: `'d42'` (default; `DataFrame.toByteArray()` output)
   * or `'parquet'` (converted client-side via the Arrow package — fails with a clear error when
   * Arrow is not installed). Ignored for other payloads — the format is inferred:
   * DataFrame → sent as d42, string → `'csv'`, object[] → `'json'`. */
  format?: 'csv' | 'd42' | 'parquet' | 'json';
}

/** Batch upload report of domain table `batch`. */
export interface DomainBatchReport {
  inserted: number;
  updated: number;
  skipped: number;
  errorCount: number;
  /** Per-row outcomes; capped server-side. */
  rows: DomainBatchRowResult[];
  /** Set when the batch failed but a per-row report is available (e.g. an allOrNothing abort). */
  error?: string;
}

/** Insert payload for domain table `insert`: row values plus an optional
 * idempotency key (for tables that declare `"idempotency": true`). */
export type DomainRowInsert<TRow> = Partial<TRow> & {idempotencyKey?: string};

/** A saved filter preset of a domain table — a small shareable entity carrying the filter
 * panel's state maps verbatim (see `DomainSavedFiltersClient`). */
export interface DomainSavedFilterInfo {
  id: string;
  name: string;
  friendlyName: string;
  /** Per-column filter states (the shape `DG.FilterGroup` saves), stored verbatim. */
  states: {[column: string]: any};
  /** The preset's author (preserved on in-place updates). */
  author?: any;
}

export interface DomainTableClientOptions {
  /** Datetime columns to materialize as dayjs on JSON reads (generated clients pass this;
   * untyped clients keep ISO strings). Dotted `'<fk>.<col>'` entries cover master-expand
   * fields. */
  datetimeColumns?: string[];
  /** Datetime columns of `'details:'` child rows, keyed by the result field (the child-table
   * name): materialized as dayjs recursively when the expand is requested. */
  detailDatetimeColumns?: {[detailField: string]: string[]};
}

export type DomainTxValues<T> = {[K in keyof T]: T[K] | `$${string}`};

/** Result type of one transaction op, keyed on its `op` discriminant — powers the
 * mapped-tuple `transaction()` signatures (per-op result types, no positional casts). */
export type DomainOpResultFor<TOp> =
  TOp extends {op: 'insert'} ? DomainInsertResult :
  TOp extends {op: 'update'} ? DomainUpdateResult :
  TOp extends {op: 'delete'} ? DomainDeleteResult : DomainOpResult;

/** Minimal structural client the builder executes against (avoids a domains.ts → dapi.ts cycle). */
export interface IDomainQueryExecutor<TRow> {
  query(spec: any): Promise<TRow[]>;
  queryDf(spec: any): Promise<any>;       // DG.DataFrame
  count(filter?: any): Promise<number>;
}

/** Bound condition node: `cond('name', '=', "O'Brien")`. The value travels in the condition
 * tree and is bound server-side — NEVER interpolated into a filter string, so any string
 * value is expressible (apostrophes included, which the smart-filter grammar cannot quote). */
export function cond<TColumn extends string = string>(
  property: DomainColumnRef<TColumn>, operator: DomainConditionOperator,
  value?: DomainFilterValue): DomainCondition<TColumn> {
  return value === undefined ? {property, operator} : {property, operator, value};
}

/** AND-combined condition tree: `and(a, b, c)` → `[a, 'and', b, 'and', c]`. */
export function and<TColumn extends string = string>(
  ...nodes: (DomainCondition<TColumn> | DomainConditionTree<TColumn>)[]): DomainConditionTree<TColumn> {
  return _joinNodes(nodes, 'and');
}

/** OR-combined condition tree: `or(a, b, c)` → `[a, 'or', b, 'or', c]`. */
export function or<TColumn extends string = string>(
  ...nodes: (DomainCondition<TColumn> | DomainConditionTree<TColumn>)[]): DomainConditionTree<TColumn> {
  return _joinNodes(nodes, 'or');
}

function _joinNodes<TColumn extends string>(
  nodes: (DomainCondition<TColumn> | DomainConditionTree<TColumn>)[],
  connector: 'and' | 'or'): DomainConditionTree<TColumn> {
  const tree: DomainConditionTree<TColumn> = [];
  for (const n of nodes) {
    if (tree.length > 0)
      tree.push(connector);
    tree.push(n as any);
  }
  return tree;
}

/** Thenable query builder (returned by no-arg `query()`): build with
 * `.where/.orderBy/.select/.expand/.top/.skip`, then `await` it (rows), or finish with
 * `.df()/.first()/.count()/.exists()`. Prefer the condition forms of `where` (and the
 * `cond`/`and`/`or` helpers) over template-built filter strings — condition values are
 * bound server-side, so any string value is safe (apostrophes included). Without `.top()`
 * the server's default limit (100) applies; page larger sets with `.top()/.skip()`.
 * Immutable-ish: `expand()`/`select()` return a re-typed builder; other methods mutate
 * and return this. */
export class DomainQueryBuilder<TRow, TColumn extends string = string,
    TExpand extends {[key: string]: {}} = {[key: string]: {}}, TResult = TRow>
    implements PromiseLike<TResult[]> {
  private _conds: DomainConditionTree<TColumn> = [];
  private _rawFilter?: string;
  private _orders: string[] = [];
  private _columns?: string[];
  private _expand: string[] = [];
  private _limit?: number;
  private _offset?: number;

  constructor(private readonly client: IDomainQueryExecutor<TRow>) {}

  /** Equality map (AND-combined), a single typed condition (3-arg or node), or a raw
   * smart-filter string escape hatch; multiple where() calls AND-combine. A raw string
   * cannot be combined with conditions (the string parses server-side) — use
   * `cond()/and()/or()` to express everything as one tree instead. */
  where(values: {[K in TColumn]?: DomainFilterValue}): this;
  where(property: DomainColumnRef<TColumn>, operator: DomainConditionOperator,
        value?: DomainFilterValue): this;
  where(filter: DomainFilter<TColumn>): this;
  where(a: any, operator?: DomainConditionOperator, value?: DomainFilterValue): this {
    if (typeof a === 'string' && operator !== undefined)
      this._addCond(cond(a, operator, value) as DomainCondition<TColumn>);
    else if (typeof a === 'string') {
      if (this._rawFilter !== undefined || this._conds.length > 0)
        throw new Error('combine string filters with conditions via cond()/and()/or()');
      this._rawFilter = a;
    }
    else if (Array.isArray(a) || (a != null && 'property' in a && 'operator' in a))
      this._addCond(a);
    else if (a != null)
      for (const k of Object.keys(a))
        this._addCond({property: k as any, operator: '=', value: a[k]});
    return this;
  }

  private _addCond(node: DomainCondition<TColumn> | DomainConditionTree<TColumn>): void {
    if (this._rawFilter !== undefined)
      throw new Error('combine string filters with conditions via cond()/and()/or()');
    if (this._conds.length > 0)
      this._conds.push('and');
    this._conds.push(node as any);
  }

  /** Appends a sort key (the `'col,!col2'` grammar; [desc] adds the `!`). */
  orderBy(column: DomainColumnRef<TColumn>, desc: boolean = false): this {
    this._orders.push(`${desc ? '!' : ''}${column}`);
    return this;
  }

  /** Narrows the projection to [columns] (system columns always ride along). */
  select<K extends TColumn & keyof TRow & string>(...columns: K[]):
      DomainQueryBuilder<TRow, TColumn, TExpand,
        Pick<TRow, K | Extract<keyof TRow, 'id' | 'version' | 'created_on' | 'updated_on' | 'author_id'>>> {
    this._columns = columns;
    return this as any;
  }

  /** Adds an expand and intersects its fields into the awaited row type
   * (`'details:'` child rows are full child rows — dayjs datetimes included). */
  expand<K extends keyof TExpand & string>(key: K):
      DomainQueryBuilder<TRow, TColumn, TExpand, TResult & TExpand[K]> {
    this._expand.push(key);
    return this as any;
  }

  /** Row cap for this query (server default 100 without it). */
  top(count: number): this {
    this._limit = count;
    return this;
  }

  /** Row offset (pair with {@link top} to page). */
  skip(count: number): this {
    this._offset = count;
    return this;
  }

  private _filter(): DomainFilter<TColumn> | undefined {
    return this._rawFilter !== undefined ? this._rawFilter
      : this._conds.length === 0 ? undefined : this._conds;
  }

  private _spec(): any {
    const spec: any = {};
    const filter = this._filter();
    if (filter !== undefined)
      spec.filter = filter;
    if (this._orders.length > 0)
      spec.sort = this._orders.join(',');
    if (this._columns != null)
      spec.columns = this._columns;
    if (this._expand.length > 0)
      spec.expand = this._expand;
    if (this._limit != null)
      spec.limit = this._limit;
    if (this._offset != null)
      spec.offset = this._offset;
    return spec;
  }

  private _run(): Promise<TResult[]> {
    return this.client.query(this._spec()) as Promise<any>;
  }

  /** The same query as a typed DataFrame (d42; `'details:'` expand is JSON-only and
   * rejected server-side — await the builder for detail arrays instead). */
  df(): Promise<any> {
    return this.client.queryDf(this._spec());
  }

  /** First matching row or null (forces `top(1)`). */
  async first(): Promise<TResult | null> {
    this._limit = 1;
    const rows = await this._run();
    return rows.length === 0 ? null : rows[0];
  }

  /** Matching-row count under the built filter (ignores top/skip/select/expand). */
  count(): Promise<number> {
    return this.client.count(this._filter());
  }

  /** Whether at least one row matches. */
  async exists(): Promise<boolean> {
    return (await this.count()) > 0;
  }

  then<TR1 = TResult[], TR2 = never>(
    onfulfilled?: ((value: TResult[]) => TR1 | PromiseLike<TR1>) | null,
    onrejected?: ((reason: any) => TR2 | PromiseLike<TR2>) | null): Promise<TR1 | TR2> {
    return this._run().then(onfulfilled, onrejected);
  }
}
