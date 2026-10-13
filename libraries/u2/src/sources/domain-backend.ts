/* The platform seam of `DomainSource` (GOAL "What it looks like"): one table handle answering
   its shape, its access, its rows and its writes, and one `saveAll` for every writer's batch as
   a single transaction. dg fills `backends.domain` over `grok.dapi.domains`; `MemoryDomainBackend`
   answers it from a schema.json in memory, so the gallery, the headless tests and an agent
   iterating on an app need no server. What both must agree on is written down in
   docs/domain-backend-contract.md. */
import type {IProperty} from '../core/property-like.js';
import type {AccessData} from '../core/access.js';
import type {DomainCondition, DomainConditionTree} from '../core/filter/model.js';
import type {DataFrameLike} from './df-bindings.js';
import type {EditState} from './edit-state.js';

/** The registry's reflection of a table (js-api `DomainTableInfo`): naming, keys, and what the
 * schema declares for search, constraints, ref filters, custom permissions and the tables that
 * reference this one. */
export interface DomainTableInfoLike {
  nameColumn: string | null;
  singularName: string;
  pluralName: string;
  /** The table's declared display name, when the registry carries one — what a view title and a
   * breadcrumb root read before falling back to the plural name. */
  friendlyName?: string;
  businessKey: string[];
  /** How an app spells a row in a URL: the business-key values, or the canonical row id as one
   * opaque segment (a key whose values may carry any delimiter). */
  rowAddress: 'id' | 'businessKey';
  /** What a query `search` runs over: the `searchable` columns, else the name column. */
  searchableColumns: string[];
  /** The grammar constraints of the schema (`{expr}` entries; SQL `check`s stay server-side). */
  constraints: {name: string, expr: string, message?: string}[];
  /** Per ref column, the grammar filter narrowing its candidates (`$param` = a sibling column). */
  refFilters: Record<string, string>;
  /** The custom permission names the schema declares — `can.<name>` / `~can_<name>`. */
  permissions: string[];
  /** Every table with a ref column pointing here; `label` is that column's caption. */
  childTables: {schema: string, table: string, fkColumn: string, label: string}[];
  /** Whether the schema declares the table a hierarchy (`hierarchy: true`): exactly one ref
   * column targets the table itself, and {@link DomainTableLike.ancestors} walks it. */
  hierarchy?: boolean;
  /** That self-referencing column, when the table is a hierarchy. */
  parentColumn?: string | null;
}

/** Which rows a query answers: the live ones (the default), the live and the soft-deleted, or
 * the deleted alone — a trash list, whose rows carry `~is_deleted` and are read-only until
 * {@link DomainTableLike.restore} brings them back. */
export type DomainDeletedMode = 'exclude' | 'include' | 'only';

/** The `DomainQuerySpec` subset a source issues: a smart-filter string or the canonical
 * condition tree, a case-insensitive `search` over the searchable columns, `'col,!col'`
 * ordering, paging, and whether the rows carry the per-row access columns (`Access.ROW_COLUMNS`). */
export interface DomainQueryLike {
  filter?: string | DomainCondition | DomainConditionTree;
  search?: string;
  sort?: string;
  columns?: string[];
  limit?: number;
  offset?: number;
  withAccess?: boolean;
  /** Default `'exclude'`; anything else projects `~is_deleted` with every row. */
  deleted?: DomainDeletedMode;
  /** Ref columns whose target's display name should ride with the rows: each projects one
   * `~caption_<column>` string cell (null where the caller may not see the target). Never a target
   * field — a caller that wants fields uses a query of its own. */
  captions?: string[];
}

/** What a read is scoped to — the part of a query that selects rows rather than shapes the page.
 * Every read takes the same object, so a caller cannot forward `filter` and forget `deleted`. */
export type DomainReadScope = Pick<DomainQueryLike, 'filter' | 'search' | 'deleted'>;

/** Whether a read selects nothing at all — the whole live table, which a change token answers
 * for. An empty condition tree is what an empty filter builder compiles to. */
export function isUnscoped(scope: DomainReadScope): boolean {
  const filter = scope.filter;
  return (filter === undefined || filter === '' || (Array.isArray(filter) && filter.length === 0)) &&
    (scope.search === undefined || scope.search === '') &&
    (scope.deleted === undefined || scope.deleted === 'exclude');
}

/** What a table can do AT ALL, independent of the caller — the backend's own answer for every
 * optional behaviour a control would otherwise guess from the table's shape. The platform
 * computes it once per handle (`DomainAccess.support`); the memory backend computes it from the
 * schema.
 *
 * The rule, once: an optional member of {@link DomainTableLike} is installed only when its flag
 * is true, and a caller checks `=== undefined` and refuses BY NAME. Nothing is guessed, and a
 * member that is there always works. */
export interface DomainSupportLike {
  /** The system columns the table physically has, in projection order — a registration declaring
   * a subset lists only what it carries. */
  systemColumns: string[];
  /** The engine accepts writes: gates {@link DomainTableLike.updateWhere} and
   * {@link DomainTableLike.batch}. */
  writes: boolean;
  /** Soft delete: `deleted: 'include' | 'only'` reads and `~is_deleted`. */
  deleted: boolean;
  /** Gates {@link DomainTableLike.restore} — and with it a `deleted` source at all. */
  restore: boolean;
  /** Gates {@link DomainTableLike.audit}. */
  audit: boolean;
  /** The table declares a hierarchy, so {@link DomainTableLike.ancestors} answers. */
  ancestors: boolean;
  /** Gates {@link DomainTableLike.probe} — the table carries `updated_on`. */
  probe: boolean;
  /** The backend keeps a change token that moves with every write, so an unscoped probe can read
   * it instead of aggregating. */
  version: boolean;
  /** The backend can tell a client that rows changed without being asked. The memory backend has
   * no subscriptions and says so. */
  watch: boolean;
  /** Gates {@link DomainTableLike.updateWhere}. */
  updateWhere: boolean;
  /** A query's `captions` project `~caption_<column>` with the rows; false where a source resolves
   * ref captions per row instead. */
  captions: boolean;
  /** {@link DomainTableLike.transaction} lands on this table — the whole write path of a source. */
  transaction: boolean;
  /** How an update is guarded against a concurrent edit: by the row `version`, by the `expected`
   * old values of the changed columns, or not at all. */
  concurrency: 'version' | 'expected' | 'none';
  /** The filter grammar the backend answers: `'basic'` has no `under`, no regex, no `!like`, no
   * datetime `!=` and no bool null test. */
  filters: 'full' | 'basic';
  /** What {@link DomainTableLike.batch} does beyond a plain insert; `validate` gates
   * {@link DomainTableLike.validate}. */
  batch: {upsert: boolean, partial: boolean, validate: boolean, skipDuplicates: boolean};
}

/** What one poll of a live source learns: how many rows match the query and when the newest of
 * them was last written — the platform's aggregate `count` + `max(updated_on)` in one request.
 * `last` is null over an empty match.
 *
 * `count: -1` means NOT COUNTED: the backend answered the whole question with a change token
 * ({@link DomainTableLike.probe} over an unscoped read), so `last` alone moves. A source only
 * ever compares the pair against its own previous poll, and re-baselines whenever the scope
 * changes, so the two shapes never meet in one comparison. */
export interface DomainProbeLike {
  count: number;
  last: string | null;
}

/** One `/transaction` op (js-api `DomainTransactionOp`): `table` is `'<table>'` in the writer's
 * schema or `'<schema>.<table>'`; an insert may name itself with `ref`, and any op's value
 * `'$<ref>'` — earlier or later in the batch — is replaced by that row's id. */
export interface DomainTransactionOpLike {
  /** A `restore` carries `id` alone and undoes a landed soft delete — the Delete grant, not a new one. */
  op: 'insert' | 'update' | 'delete' | 'restore';
  table: string;
  ref?: string;
  values?: Record<string, unknown>;
  id?: string;
  expectedVersion?: number;
  /** The old-value guard of an update where `support.concurrency` is `'expected'`: the columns'
   * values as last read; never together with `expectedVersion`. */
  expected?: Record<string, unknown>;
}

export interface DomainTransactionResultLike {
  id?: string;
  version?: number;
  /** The system columns the write stamped. The platform's `/transaction` answers `{id, version}`
   * alone and the js-api editor re-reads the row for the rest; the memory backend has the stamped
   * row at hand and answers with it, so both writers leave the same frame behind. */
  created_on?: string;
  updated_on?: string;
  author_id?: string;
}

/** What a bulk upload does with rows the table already has (js-api `DomainBatchOptions`). */
export interface DomainBatchOptionsLike {
  /** `'insert'` (the default) or `'upsert'` — merge by the table's business key. */
  mode?: 'insert' | 'upsert';
  /** Abort the whole upload on any row error (default true). */
  allOrNothing?: boolean;
  /** Report business-key duplicates as errors instead of skipping them. */
  errorOnDuplicate?: boolean;
}

/** What a bulk upload answers (js-api `DomainBatchReport`): the counts, one line per row, and
 * `error` where the upload failed but a per-row report survived it. A storage that cannot tell
 * an insert from an update answers an upsert as `merged` (per row `status: 'merged'`). */
/** The totals are absent on a warehouse refusal (a constraint the connector reports), which
 * carries `error` and the failing rows alone. */
export interface DomainBatchReportLike {
  inserted?: number;
  updated?: number;
  merged?: number;
  skipped?: number;
  errorCount?: number;
  rows: {index: number, id: string | null, status: string, existingId?: string,
    errors?: {column?: string, code?: string, message: string}[]}[];
  error?: string;
}

/** One row of a dry run: what the batch WOULD do with it. */
export interface DomainBatchValidationRowLike {
  index: number;
  predicted: 'insert' | 'update' | 'skip' | 'error';
  /** The row the prediction is about, for `update` and `skip` — a predicted insert has no id. */
  existingId?: string;
  errors?: {column?: string, code?: string, message: string}[];
}

/** What a dry run answers (js-api `DomainBatchValidation`): the counts and one line per row, and
 * no ids or statuses — nothing was written, and a verdict is about the table AS IT IS NOW. */
export interface DomainBatchValidationLike {
  validateOnly: true;
  rowCount: number;
  willInsert: number;
  willUpdate: number;
  willSkip: number;
  errorCount: number;
  rows: DomainBatchValidationRowLike[];
}

/** One line of a row's history (js-api `DomainAuditEntry`): what an op did to it, under which
 * transaction, by whom. */
export interface AuditEntryLike {
  id: string;
  tx_id: string;
  op: string;
  actor_id: string | null;
  ts: string;
  before: Record<string, unknown> | null;
  after: Record<string, unknown> | null;
}

/** A frame the backend fetched together with its single writer: the frame is the collection,
 * `edit` tracks every change made to it, `append` adds the next page into the same frame.
 * Replaced whole on a re-query, so `dispose` is where the writer detaches (STATE-CONTRACT H8). */
export interface DomainFrameLike {
  df: DataFrameLike;
  edit: EditState;
  /** Answers the number of rows appended. */
  append(spec: DomainQueryLike): Promise<number>;
  dispose(): void;
}

export interface DomainTableLike {
  /** `'<schema>.<table>'`. */
  address: string;
  properties: IProperty[];
  info: DomainTableInfoLike;
  /** What the table can do at all — declared, never guessed; see {@link DomainSupportLike}. */
  readonly support: DomainSupportLike;
  access(): Promise<AccessData>;
  query(spec: DomainQueryLike): Promise<Record<string, unknown>[]>;
  /** The total under the same scope a query takes. */
  count(scope?: DomainReadScope): Promise<number>;
  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]>;
  /** The rows as a frame with the writer attached — the one collection a source holds. */
  frame(spec: DomainQueryLike): Promise<DomainFrameLike>;
  /** A row's history, oldest first — the platform's `client.audit`. */
  audit?(id: string): Promise<AuditEntryLike[]>;
  /** Brings a soft-deleted row back (`~is_deleted` off, the `'undelete'` audit op); refused with
   * `validation` where the row refers to a deleted parent. A backend that does not declare it
   * cannot answer a `deleted` query either — `DomainSource` refuses one over it. */
  restore?(id: string): Promise<void>;
  /** Writes `values` into every live row the filter matches that the caller may edit, as ONE
   * transaction (the platform's `POST …/{table}/update`): the filter is required, a column the
   * caller may not write is refused before anything is written, and at most `limit` rows are
   * touched — `hasMore` says the filter matched more than that. */
  updateWhere?(filter: DomainQueryLike['filter'], values: Record<string, unknown>,
    options?: {limit?: number}): Promise<{updated: number, hasMore: boolean}>;
  /** The row's ancestors along the table's `parentColumn`, ROOT FIRST and without the row
   * itself (the platform's `GET …/{id}/path`). Only a hierarchy table answers it; the chain
   * stops at the first ancestor the caller cannot see, at a cycle, and at depth 64. */
  ancestors?(id: string): Promise<{id: string, name: string}[]>;
  /** What the collection looks like on the server right now, under the same scope a query takes —
   * ONE request, never a page of rows: a `live` source polls it and refreshes when the pair
   * moved. A backend that does not declare it is not polled, so a `live` source over one is
   * exactly as live as the backend can be: not at all. */
  probe?(scope?: DomainReadScope): Promise<DomainProbeLike>;
  /** Uploads whole rows in one call (the platform's `POST …/{table}/batch`), which is where the
   * upsert merge, the duplicate rules and the per-row report live. A backend that does not
   * declare it cannot be imported into — `domains.import` refuses by name before the wizard
   * opens, rather than standing a transaction in for the real thing. */
  batch?(rows: Record<string, unknown>[], options?: DomainBatchOptionsLike): Promise<DomainBatchReportLike>;
  /** What {@link batch} WOULD do with these rows, without doing it: every check the commit runs —
   * coercions, the schema's rules, business-key duplicates inside the batch and against the live
   * rows, FK existence — inside a transaction that is rolled back. Installed wherever `batch` is,
   * so an import preview is the backend's verdict instead of a second implementation of it. */
  validate?(rows: Record<string, unknown>[], options?: DomainBatchOptionsLike):
    Promise<DomainBatchValidationLike>;
}

export interface DomainBackend {
  table(address: string): Promise<DomainTableLike>;
  /** The backend's batched display-name resolver for ref ids — the platform registry's, which
   * coalesces every caller's ids into one narrow fetch per table and caches them. Every requested
   * id is a key; an id the caller may not see, or that no row answers, maps to null. A backend
   * without one is read row by row. */
  resolveNames?(table: string, ids: readonly string[]): Promise<Record<string, string | null>>;
  /** Every writer's pending batch as ONE transaction; resolves to whether it landed. The memory
   * backend throws its `DomainBackendError`; the platform answers false after its own dialogs
   * (the js-api `DomainSession` owns the conflict and validation loop). */
  saveAll(edits: EditState[]): Promise<boolean>;
}

/** What a backend throws: the server's error family by `code` (`'version-conflict'`,
 * `'validation'`, `'not-found'`, `'bad-ref'`), the message for the user. */
export class DomainBackendError extends Error {
  constructor(readonly code: string, message: string) {
    super(message);
    this.name = 'DomainBackendError';
  }
}
