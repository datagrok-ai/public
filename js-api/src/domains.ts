/**
 * Typed surface for entity-mapped domain tables (`grok.dapi.domains`): condition-tree
 * filter types, typed results, the DomainError class family, predicate helpers, the fluent
 * query builder, and optimistic-concurrency helpers.
 *
 * Two rules hold everywhere: condition VALUES are bound server-side — never interpolated
 * into filter strings, so any string value is expressible (apostrophes included, which the
 * smart-filter string grammar cannot quote); and datetime columns materialize as dayjs on
 * generated clients (`datetimeColumns`) while untyped clients keep ISO strings.
 *
 * @remarks BREAKING (codegen v2, GROK-20602): `grok api`-generated clients type datetime
 * columns as `Dayjs` and thread `<Table>Column`/`<Table>Expand` generics; regenerating
 * db.ts changes its surface. The generated `dapi2.domains` namespace is removed
 * (GROK-20601) — this module and `DomainTableClient` are the API.
 */
import type {Dayjs} from 'dayjs';
import type {DataFrame} from './dataframe';

// The runtime DG namespace, declared instead of imported: dg.ts re-exports this module,
// so importing it back would close a module cycle (same seam as IDomainQueryExecutor).
declare let DG: any;

/** The system columns every domain table carries (always projected on reads). */
export type DomainSystemColumn = 'id' | 'version' | 'created_on' | 'updated_on' | 'author_id';

/** Column reference: a declared column of the table, a system column, or a dotted FK path
 * ('project_id.name', up to 3 forward hops). Used by the condition, sort, and aggregate
 * column slots — those accept system columns even when a hand-written [TColumn] union omits
 * them (generated `<Table>Column` unions include them). Other slots (equality-map `where`,
 * `columns`, `select()`, `fetchFields()`, column security) stay TColumn-typed. */
export type DomainColumnRef<TColumn extends string = string> =
  TColumn | DomainSystemColumn | `${string}.${string}`;

/** Operators of the canonical condition tree (server: filter_compiler.dart:149-153 + fuzzy). */
export type DomainConditionOperator =
  '=' | '!=' | '>' | '>=' | '<' | '<=' | 'like' | 'not like' | '~*' | '!~*' | 'is' | 'is not' | 'fuzzy';

/** Value of one condition: scalar, list (`= ANY` / `!= ALL`), null (IS [NOT] NULL),
 * or dayjs (sent as ISO-8601). Always bound server-side. */
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
/** The canonical condition tree: conditions and nested trees joined by connector strings. */
export type DomainConditionTree<TColumn extends string = string> =
  DomainConditionNode<TColumn>[];

/** Every filter slot accepts the smart-filter string (human shorthand), a single condition,
 * or the canonical condition tree — uniformly across query, queryDf, aggregate, facets,
 * count, exists, and deleteWhere (server: DomainFilterCompiler.toTree). */
export type DomainFilter<TColumn extends string = string> =
  string | DomainCondition<TColumn> | DomainConditionTree<TColumn>;

/** Runtime list of the system columns (type-level counterpart: {@link DomainSystemColumn}). */
export const DOMAIN_SYSTEM_COLUMNS = ['id', 'version', 'created_on', 'updated_on', 'author_id'] as const;

/** Splits the `'<schema>.<table>'` address every domain client and UI component
 * takes, throwing on a malformed one — the single spelling of that contract. */
export function splitDomainTable(name: string): [string, string] {
  const dot = name == null ? -1 : name.indexOf('.');
  if (dot < 1 || dot === name.length - 1)
    throw new Error(`Domain table name must be '<schema>.<table>', got '${name}'`);
  return [name.substring(0, dot), name.substring(dot + 1)];
}

/** One per-column failure of a validation-failed write. */
export interface DomainColumnError { column: string; code: string; message: string; }
/** Per-row failure list of a validation-failed write ({@link DomainValidationError.rows}). */
export interface DomainRowErrors { index: number; id?: string | null; errors: DomainColumnError[]; }

/** Result of one inserted row. */
export interface DomainInsertResult {
  id: string;
  created: boolean;
  version?: number;
  /** Set when the row was NOT created: 'duplicate' (business-key match) or 'idempotent-replay'. */
  status?: 'duplicate' | 'idempotent-replay';
  existingId?: string;
}
/** Result of one updated row: `version` is the new (incremented) row version. */
export interface DomainUpdateResult { id: string; version: number; }
/** Result of one deleted row. */
export interface DomainDeleteResult { id: string; deleted: true; }
/** Result of one transaction op (see {@link DomainOpResultFor} for the typed-tuple form). */
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

/** Grantable permission on a domain registry entity (table, schema, or column schema). */
export type DomainPermission = 'View' | 'Edit' | 'Delete' | 'Share';

/** One direct permission row on a domain registry entity (see `DomainTableClient.grants`). */
export interface DomainGrant {
  group: {id: string; friendlyName: string; personal: boolean};
  permission: DomainPermission;
}

/** Effective capabilities of the CURRENT user on one domain table
 * (see `DomainTableClient.capabilities`). Composed client-side from server-truth
 * permission probes on the FINAL securing entity through the master delegate chain,
 * plus the writable-column mirror of column security.
 *
 * These are TABLE-level affordances, and they drift toward denial: every flag is
 * derived from grants on the securing entity, so grants that reach individual ROWS
 * (a master row of a master-mode table, a promoted row of a row-mode table) are not
 * counted. A false flag therefore means "no table-wide right" — the caller may still
 * succeed on a particular row, and per-row truth comes from
 * {@link DomainRow.permissions}. A true flag mirrors a predicate the server also
 * enforces, so gating UI on it avoids the common 403s (but never assume it removes
 * them: grants can change between the probe and the write, and column-level
 * restrictions are checked per value). */
export interface DomainTableCapabilities {
  /** View grant on the securing entity. Row-mode tables may still expose
   * individually granted rows when false. */
  canView: boolean;
  /** Edit grant on the securing entity — the server's insert predicate
   * (`_checkInsertAllowed`). False negative on master-mode tables where the caller
   * holds the grant on an individual MASTER ROW rather than the master table. */
  canInsert: boolean;
  /** {@link canInsert} AND at least one column is writable for the caller (the
   * built-in grid-editability rule) AND the table grant actually reaches rows —
   * for non-table security modes that needs the securing table's rows to default
   * to table visibility, otherwise access is per-row. Same false-negative shape as
   * {@link canInsert}. */
  canEdit: boolean;
  /** Delete grant on the securing entity, under the same reaches-rows rule as
   * {@link canEdit}; same false-negative shape. */
  canDelete: boolean;
  canShareTable: boolean;
  /** Column names the caller may write (Edit on an owning property schema),
   * in declared column order. */
  writableColumns: string[];
  securityMode: 'table' | 'master' | 'row';
  /** Whether writes leave an in-transaction audit trail. */
  audit: boolean;
  hasBusinessKey: boolean;
}

/** FK-inverted reference to a child (detail) table — drives detail-table links
 * and entity-view child tabs (see `DomainRegistryClient.tableInfo`). */
export interface DomainChildTableRef {
  schema: string;
  table: string;
  /** The child table's FK column referencing this table. */
  fkColumn: string;
  /** Friendly label of the FK column ('issue_id' → 'Issue') — disambiguates
   * children referencing the same parent through several columns. */
  label: string;
}

/** Registry metadata of one domain table (see `grok.dapi.domains.registry`). */
export interface DomainTableInfo {
  /** Primary display-name column; null when the table declares none. */
  nameColumn: string | null;
  /** Natural-key columns; empty when the table declares none. */
  businessKey: string[];
  /** Effective singular display name (declared, else derived from the table name). */
  singularName: string;
  /** Effective plural display name (declared, else derived from the table name). */
  pluralName: string;
  securityMode: 'table' | 'master' | 'row';
  audit: boolean;
  childTables: DomainChildTableRef[];
}

/** Per-row outcome inside a {@link DomainBatchReport}. */
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

/** Parameter values of the `DomainQuery` function — the serializable representation of
 * what a domain view is showing (filter elements, ordering, cap), as
 * {@link DomainView.query} reports it. Feeding it back to the function
 * ({@link DomainQuery.run}) reproduces the subset as a DataFrame, with the creation
 * script recorded.
 *
 * This is the ONLY parameter shape of the function — {@link DomainQuery} is its
 * mutable, URL-serializable counterpart. Grammar strings, not condition trees: each
 * element is one `DomainQuery` filter expression (`'status = "open"'`) with the
 * per-element JSON escape hatch. Use {@link DomainQuerySpec} for the REST surface
 * instead.
 *
 * Aggregate mode (non-empty `aggregations`/`groupBy`) REJECTS `columns` and `offset`
 * with a 400 — they are select-mode parameters, not ignored ones. */
export interface DomainQueryParams {
  schema: string;
  table: string;
  /** Projection; omit for all viewable columns (select mode only). */
  columns?: string[];
  /** Filter elements, AND-combined. OMITTED, not `[]`, when there are none —
   * these are the function's parameter values, and it takes no empty lists
   * (`DomainView.query` on an unfiltered view has no `filters` key). */
  filters?: string[];
  /** Master-FK expands (`'<fk_column>'`); `'details:'` child arrays are not supported. */
  joins?: string[];
  /** Measures (`'count'`, `'avg(amount) as avg_amount'`) — non-empty means aggregate mode. */
  aggregations?: string[];
  /** Grouping columns — non-empty means aggregate mode. */
  groupBy?: string[];
  /** Sort elements; a `'!'` prefix means descending. */
  orderBy?: string[];
  limit?: number;
  /** Row offset (select mode only). */
  offset?: number;
}

/** Query options for domain table `query`/`queryDf`. */
export interface DomainQuerySpec<TColumn extends string = string, TExpandKey extends string = string> {
  /** Smart-filter string (same grammar as entity search, e.g. `barcode starts "P-1"`),
   * a single condition, or the canonical condition tree — values are bound server-side.
   * Row caps: JSON `query` 10k, d42 `queryDf` 10M. */
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
/** One result row of `aggregate` — keys are the group columns and measure aliases. */
export type DomainAggregateRow<TKeys extends string = string> =
  {[K in TKeys]: number | string | boolean | null};

/** The facet kinds of {@link DomainFacetSpec}. */
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
/** Result of a `'categories'` facet; `hasMore` set when the category cap was hit. */
export interface DomainFacetCategoriesResult { categories: DomainFacetCategory[]; hasMore?: boolean; }
/** Result of a `'histogram'` facet: `buckets` respect the filter, `totalBuckets` ignore it. */
export interface DomainFacetHistogramResult {
  min: number | string | null; max: number | string | null;
  buckets: number[]; totalBuckets: number[]; nulls: number;
}
/** Result of a `'minMax'` facet (row predicate only — the stable-axis rule ignores the filter). */
export interface DomainFacetMinMaxResult { min: number | string | null; max: number | string | null; }
/** Result of a `'count'` facet: the filtered row count. */
export interface DomainFacetCountResult { count: number; }
/** Result of a `'plan'` facet: per-column profile for choosing filter-control types. */
export interface DomainFacetPlanResult {
  columns: {name: string; distinct: number; min?: number | string; max?: number | string}[];
}
/** Maps a facet `kind` to its result type (keys the typed `facets()` response). */
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

/** Options of `DomainsDataSource.table` (generated clients pass these). */
export interface DomainTableClientOptions {
  /** Datetime columns to materialize as dayjs on JSON reads (generated clients pass this;
   * untyped clients keep ISO strings). Dotted `'<fk>.<col>'` entries cover master-expand
   * fields. */
  datetimeColumns?: string[];
  /** Datetime columns of `'details:'` child rows, keyed by the result field (the child-table
   * name): materialized as dayjs recursively when the expand is requested. */
  detailDatetimeColumns?: {[detailField: string]: string[]};
}

/** Typed transaction-op values: each column also accepts a `'$<ref>'` back-reference. */
export type DomainTxValues<T> = {[K in keyof T]: T[K] | `$${string}`};

/** Result type of one transaction op, keyed on its `op` discriminant — powers the
 * mapped-tuple `transaction()` signatures (per-op result types, no positional casts). */
export type DomainOpResultFor<TOp> =
  TOp extends {op: 'insert'} ? DomainInsertResult :
  TOp extends {op: 'update'} ? DomainUpdateResult :
  TOp extends {op: 'delete'} ? DomainDeleteResult : DomainOpResult;

/** Runs [action], retrying on DomainVersionConflictError (for transaction-based
 * read-modify-write flows — put the fresh read INSIDE [action]); rethrows anything else
 * and the final conflict. `maxRetries` counts retries AFTER the initial attempt
 * (default 5 retries = up to 6 attempts). Retries run immediately, without backoff —
 * conflicts resolve by re-reading, not waiting; add your own delay for high-contention
 * hot rows. */
export async function retryOnVersionConflict<T>(
    action: () => Promise<T>, options?: {maxRetries?: number}): Promise<T> {
  const maxRetries = options?.maxRetries ?? 5;
  for (let attempt = 0; ; attempt++) {
    try {
      return await action();
    } catch (e) {
      if (attempt >= maxRetries || !(e instanceof DomainVersionConflictError))
        throw e;
    }
  }
}

/** Minimal structural client the builder executes against (avoids a domains.ts → dapi.ts cycle). */
export interface IDomainQueryExecutor<TRow> {
  query(spec: any): Promise<TRow[]>;
  queryDf(spec: any): Promise<any>;       // DG.DataFrame
  count(filter?: any): Promise<number>;
  /** Table address, needed by {@link DomainQueryBuilder.toQuery} (`DomainTableClient` carries it). */
  readonly schema?: string;
  readonly table?: string;
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
 * and return this. Invariants: the builder is PromiseLike only — there is no
 * `.catch()`/`.finally()`, use `try { await b } catch`; EVERY `await` (or terminal)
 * re-executes the query — two awaits are two round trips, cache the rows instead. */
export class DomainQueryBuilder<TRow, TColumn extends string = string,
    TExpand extends {[key: string]: {}} = {[key: string]: {}}, TResult = TRow, TDf = any>
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
        throw new Error('cannot combine a string filter with conditions — express everything as one tree via cond()/and()/or()');
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
      throw new Error('cannot combine a string filter with conditions — express everything as one tree via cond()/and()/or()');
    if (this._conds.length > 0)
      this._conds.push('and');
    this._conds.push(node as any);
  }

  /** Appends a sort key (the `'col,!col2'` grammar; [desc] adds the `!`). */
  orderBy(column: DomainColumnRef<TColumn>, desc: boolean = false): this {
    this._orders.push(`${desc ? '!' : ''}${column}`);
    return this;
  }

  /** Narrows the projection to [columns] (system columns always ride along).
   * NB: select() re-types the result from TRow — call it BEFORE expand(): the
   * select-then-expand order composes, while expand-then-select keeps the expanded
   * fields at runtime but drops them from the compile-time type. */
  select<K extends TColumn & keyof TRow & string>(...columns: K[]):
      DomainQueryBuilder<TRow, TColumn, TExpand,
        Pick<TRow, K | Extract<keyof TRow, DomainSystemColumn>>, TDf> {
    this._columns = columns;
    return this as any;
  }

  /** Adds an expand and intersects its fields into the awaited row type
   * (`'details:'` child rows are full child rows — dayjs datetimes included). */
  expand<K extends keyof TExpand & string>(key: K):
      DomainQueryBuilder<TRow, TColumn, TExpand, TResult & TExpand[K], TDf> {
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
  df(): Promise<TDf> {
    return this.client.queryDf(this._spec());
  }

  /** First matching row or null. NB: permanently overwrites the builder's limit with 1 —
   * a later `await` of the same builder returns at most one row. */
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

  /** This builder's accumulated state as a serializable {@link DomainQuery} (the
   * URL / deep-link / recorded-run form of the same query). Conditions become filter
   * elements: top-level AND conjuncts are emitted one per element (so a URL can bind
   * `filters[0]` alone), anything else travels as one JSON element.
   *
   * The row cap travels, but its DEFAULT does not: without `.top()`, awaiting the
   * builder takes the server default of 100 rows while `toQuery().toSpec()` falls back
   * to {@link DOMAIN_QUERY_ROW_LIMIT} — the same builder, two row counts. Pin one with
   * `.top()` when the two forms must agree. */
  toQuery(): DomainQuery {
    if (this.client.schema == null || this.client.table == null)
      throw new Error('the query builder has no table address — construct the DomainQuery explicitly');
    return new DomainQuery({
      schema: this.client.schema, table: this.client.table,
      filters: this._rawFilter !== undefined ? [this._rawFilter] : _treeToFilterElements(this._conds),
      columns: this._columns, joins: this._expand, orderBy: this._orders,
      limit: this._limit, offset: this._offset,
    });
  }

  then<TR1 = TResult[], TR2 = never>(
    onfulfilled?: ((value: TResult[]) => TR1 | PromiseLike<TR1>) | null,
    onrejected?: ((reason: any) => TR2 | PromiseLike<TR2>) | null): Promise<TR1 | TR2> {
    return this._run().then(onfulfilled, onrejected);
  }
}

/** The `DomainQuery` function's client-side row ceiling (the server clamps to it):
 * the limit {@link DomainQuery.toSpec} writes when the query declares none, so a spec
 * never silently inherits the server's default of 100. */
export const DOMAIN_QUERY_ROW_LIMIT = 10000000;

/** List-typed `DomainQuery` parameters, in the order {@link DomainQuery.toUrlParams} emits them. */
const _DOMAIN_QUERY_LISTS = ['columns', 'filters', 'joins', 'aggregations', 'groupBy', 'orderBy'] as const;
type _DomainQueryList = (typeof _DOMAIN_QUERY_LISTS)[number];

function _isListParam(name: string): name is _DomainQueryList {
  return (_DOMAIN_QUERY_LISTS as readonly string[]).includes(name);
}

/** Copy of a list parameter; an empty or absent list normalizes to undefined, so
 * every DomainQuery has exactly ONE representation of "no elements". */
function _copyList(list?: string[] | null): string[] | undefined {
  return list == null || list.length === 0 ? undefined : list.slice();
}

/** Condition tree → `filters` elements: top-level AND conjuncts become one element each
 * (a URL can then bind `filters[0]` alone); anything with an 'or' at the top travels as
 * a single `'['`-prefixed JSON sub-group, which the function splices verbatim. */
function _treeToFilterElements(tree: DomainConditionTree): string[] | undefined {
  if (tree == null || tree.length === 0)
    return undefined;
  if (tree.some((n) => n === 'or'))
    return [JSON.stringify(tree)];
  return tree.filter((n) => typeof n !== 'string').map((n) => JSON.stringify(n));
}

/** One `filters` element as a condition node, or undefined when it is a smart-filter
 * grammar string (which only the function parses — client-side, via the Dart parser). */
function _decodeFilterElement(element: string): any {
  const s = `${element}`.trim();
  if (s.startsWith('{')) {
    let node: any;
    try {
      node = JSON.parse(s);
    } catch (_) {
      throw new Error(`Cannot parse filter "${element}": invalid JSON`);
    }
    if (node == null || typeof node !== 'object')
      throw new Error(`Cannot parse filter "${element}": expected a condition node`);
    return node;
  }
  // A '['-element that decodes as JSON is a condition sub-group; one that does not
  // is grammar ('[bracketed column] = ...').
  if (s.startsWith('[')) {
    let decoded: any;
    try {
      decoded = JSON.parse(s);
    } catch (_) { /* falls through to the grammar */ }
    if (Array.isArray(decoded)) {
      // Same minimal shape check the function applies (domain_query_func.dart): members
      // are conditions, nested groups, or 'and'/'or' — a bare-string member would splice
      // server-side nonsense.
      if (decoded.length === 0 || !decoded.every((e: any) =>
        (e != null && typeof e === 'object') || e === 'and' || e === 'or'))
        throw new Error(`Cannot parse filter "${element}": expected a condition sub-group`);
      return decoded;
    }
  }
  return undefined;
}

/** `limit`/`offset` from a URL: a negative value is rejected rather than passed on —
 * the server clamps it to 0, which would silently return nothing. */
function _parseUrlInt(key: string, value: string): number {
  const s = `${value}`.trim();
  if (!/^\d+$/.test(s))
    throw new Error(`Malformed URL parameter '${key}': expected a non-negative integer, got '${value}'`);
  return parseInt(s, 10);
}

/**
 * What the user is looking at, as ONE serializable object: the parameters of the
 * platform's `DomainQuery` function (filters, joins, aggregations, groupBy, orderBy,
 * projection, limit, offset). URL routing, deep links, saved filters, "open in Table
 * View" and data-synced dashboards all serialize exactly this — there is no parallel
 * query-state vocabulary.
 *
 * ```ts
 * const q = new DG.DomainQuery({schema: 'grit', table: 'issue',
 *   filters: ['status = "open"'], orderBy: ['!created_on'], limit: 100});
 * const df = await q.run();                                  // recorded run
 * const url = new URLSearchParams(q.toUrlParams()).toString(); // filters[0]=...&orderBy[0]=...
 * ```
 *
 * **UI-only state stays out of it.** Search text, view mode, the current entity ride
 * separate reserved URL parameters (`view=`, `entity=`) that {@link fromUrlParams}
 * ignores; if an app needs to carry them together, wrap a DomainQuery in an envelope
 * object — never widen this class.
 */
export class DomainQuery {
  schema: string;
  table: string;
  /** Projection; omit for all viewable columns (select mode only). */
  columns?: string[];
  /** Filter elements, AND-combined: smart-filter grammar strings (`'status = "open"'`),
   * or the per-element JSON escape hatch (a `'{'`-prefixed condition node, a
   * `'['`-prefixed condition sub-group). Values inside a JSON node are bound
   * server-side; a grammar string cannot quote apostrophes, so prefer nodes for
   * arbitrary user values. */
  filters?: string[];
  /** Master-FK expands (`'<fk_column>'`); `'details:'` child arrays are rejected. */
  joins?: string[];
  /** Measures (`'count'`, `'avg(amount) as avg_amount'`) — non-empty means aggregate mode. */
  aggregations?: string[];
  /** Grouping columns — non-empty means aggregate mode. */
  groupBy?: string[];
  /** Sort elements; a `'!'` prefix means descending. */
  orderBy?: string[];
  /** Row cap; {@link toSpec} falls back to {@link DOMAIN_QUERY_ROW_LIMIT} when unset. */
  limit?: number;
  /** Row offset (select mode only). */
  offset?: number;

  /** Empty and absent lists are the same thing: they normalize to undefined, so
   * {@link toParams} / {@link toUrlParams} round-trip to an identical object. */
  constructor(params: DomainQueryParams) {
    if (params == null || params.schema == null || params.table == null)
      throw new Error("DomainQuery needs a 'schema' and a 'table'");
    this.schema = params.schema;
    this.table = params.table;
    for (const name of _DOMAIN_QUERY_LISTS)
      this[name] = _copyList(params[name]);
    if (params.limit != null)
      this.limit = params.limit;
    if (params.offset != null)
      this.offset = params.offset;
  }

  /** The query behind a view's current state: `DG.DomainQuery.fromParams(domainView.query)`. */
  static fromParams(params: DomainQueryParams): DomainQuery {
    return new DomainQuery(params);
  }

  /** The state a {@link DomainQueryBuilder} accumulated (see its `toQuery`). */
  static fromBuilder(builder: DomainQueryBuilder<any>): DomainQuery {
    return builder.toQuery();
  }

  /** The `DomainQuery` function's parameter values — the shape {@link DomainView.query}
   * reports and {@link run} passes to the function. */
  toParams(): DomainQueryParams {
    const params: DomainQueryParams = {schema: this.schema, table: this.table};
    for (const name of _DOMAIN_QUERY_LISTS) {
      const list = _copyList(this[name]);
      if (list !== undefined)
        params[name] = list;
    }
    if (this.limit != null)
      params.limit = this.limit;
    if (this.offset != null)
      params.offset = this.offset;
    return params;
  }

  /** Aggregate mode: `aggregations` or `groupBy` carries elements. `columns` and
   * `offset` are then REJECTED, not ignored — {@link toSpec} throws and {@link run}
   * gets a 400 from the function. */
  get isAggregate(): boolean {
    return (this.aggregations?.length ?? 0) > 0 || (this.groupBy?.length ?? 0) > 0;
  }

  /** Runs the query through the platform's `DomainQuery` function — the reproducible
   * path: the resulting DataFrame carries a creation script, so it refreshes from the
   * Source pane, survives a saved project as a data-synced dashboard, and takes URL
   * parameters (`filters[0]`, ...). Being a normal function run, its result also goes
   * through the platform's default handling: the frame is added to the workspace and
   * OPENED IN A TABLE VIEW, which becomes the current view. Per-caller row and column
   * security applies on every run.
   *
   * The creation script is recorded only when the user has data history enabled
   * (`grok.shell.settings.dataHistory`, on by default) — the query itself runs either
   * way, but the frame comes back without `.script`/`.history` tags when it is off.
   *
   * Failures arrive as the function's own error (raw Dart message text), NOT as the
   * typed {@link DomainError} family: the function boundary carries no error envelope.
   * Where typed failures matter, run {@link toSpec} through `queryDf` instead.
   *
   * For a silent read (no history, no workspace entry, no view) pass {@link toSpec} to
   * `grok.dapi.domains.table('<schema>.<table>').queryDf(...)` instead. */
  async run(): Promise<DataFrame> {
    const func = DG.Func.byName('DomainQuery');
    if (func == null)
      throw new Error("the 'DomainQuery' function is not registered on this client");
    const call = await func.prepare(this.toParams()).call(false, undefined, {processed: false});
    return call.getOutputParamValue();
  }

  /** The same query as a REST spec for `queryDf`/`query` — the same rows as {@link run},
   * without the recording (the limit is always explicit, so the result never depends on
   * the server default).
   *
   * Throws for an aggregate query (use {@link run}, or `aggregateDf` with a
   * {@link DomainAggregateSpec}), and for several smart-filter grammar elements at
   * once: those are parsed by the function itself, and a bare string inside a REST
   * condition tree means a connector, not a filter. */
  toSpec(): DomainQuerySpec {
    if (this.isAggregate)
      throw new Error('an aggregate DomainQuery has no queryDf spec — use run(), or ' +
        'grok.dapi.domains.table(...).aggregateDf() with a DomainAggregateSpec');
    const spec: DomainQuerySpec = {limit: this.limit ?? DOMAIN_QUERY_ROW_LIMIT, offset: this.offset ?? 0};
    const filter = this._toFilter();
    if (filter !== undefined)
      spec.filter = filter;
    if (this.orderBy != null)
      spec.sort = this.orderBy.join(',');
    if (this.joins != null)
      spec.expand = this.joins;
    if (this.columns != null)
      spec.columns = this.columns;
    return spec;
  }

  private _toFilter(): DomainFilter | undefined {
    const elements = (this.filters ?? []).filter((e) => `${e}`.trim() !== '');
    if (elements.length === 0)
      return undefined;
    const nodes: DomainConditionTree = [];
    for (const element of elements) {
      const node = _decodeFilterElement(element);
      if (node === undefined) {
        if (elements.length > 1)
          throw new Error(`toSpec() cannot combine ${elements.length} filter elements into one REST ` +
            `filter: "${element}" is a smart-filter string, and only the DomainQuery function parses ` +
            'those. Use run(), or express the filters as condition nodes (cond()/and()/or()).');
        return `${element}`.trim();
      }
      if (nodes.length > 0)
        nodes.push('and');
      nodes.push(node);
    }
    return nodes;
  }

  /** URL parameters for a deep link, binding one list element per key
   * (`filters[0]=status %3D "open"`) — the platform's list-element binding scheme, so a
   * recorded run's `filters[0]` can be substituted from the URL. The table address
   * (`schema`/`table`) is NOT emitted: it addresses the view, and {@link fromUrlParams}
   * takes it back explicitly. Values are raw — encode them (`URLSearchParams`). */
  toUrlParams(): {[key: string]: string} {
    const params: {[key: string]: string} = {};
    for (const name of _DOMAIN_QUERY_LISTS) {
      const list = this[name];
      if (list != null)
        for (let i = 0; i < list.length; i++)
          params[`${name}[${i}]`] = list[i];
    }
    if (this.limit != null)
      params['limit'] = `${this.limit}`;
    if (this.offset != null)
      params['offset'] = `${this.offset}`;
    return params;
  }

  /** Rebuilds the query from {@link toUrlParams} output (lossless round trip). Keys that
   * are not query parameters — the reserved `view=` / `entity=` UI state among them —
   * are ignored; a malformed element index (`filters[x]`) or a non-integer
   * `limit`/`offset` throws instead of degrading to NaN. Element indices are read in
   * ascending order and gaps are closed. */
  static fromUrlParams(schema: string, table: string, params: {[key: string]: string}): DomainQuery {
    const query = new DomainQuery({schema, table});
    const lists: {[name: string]: {index: number, value: string}[]} = {};
    for (const key of Object.keys(params ?? {})) {
      const value = params[key];
      const bracket = key.indexOf('[');
      const name = bracket < 0 ? key : key.substring(0, bracket);
      if (bracket < 0) {
        if (name === 'limit')
          query.limit = _parseUrlInt(key, value);
        else if (name === 'offset')
          query.offset = _parseUrlInt(key, value);
        else if (_isListParam(name))
          throw new Error(`Malformed URL parameter '${key}': '${name}' is a list — ` +
            `address its elements as '${name}[0]'`);
        continue;
      }
      if (!_isListParam(name))
        continue;
      const index = /^\[(\d+)]$/.exec(key.substring(bracket));
      if (index == null)
        throw new Error(`Malformed URL parameter '${key}': expected '${name}[<index>]'`);
      if (lists[name] == null)
        lists[name] = [];
      lists[name].push({index: parseInt(index[1], 10), value: value});
    }
    for (const name of Object.keys(lists))
      query[name as _DomainQueryList] = lists[name].sort((a, b) => a.index - b.index).map((e) => e.value);
    return query;
  }
}
