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
   * untyped clients keep ISO strings). */
  datetimeColumns?: string[];
}

export type DomainTxValues<T> = {[K in keyof T]: T[K] | `$${string}`};
