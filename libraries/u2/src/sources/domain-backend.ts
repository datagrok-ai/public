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
}

/** One `/transaction` op (js-api `DomainTransactionOp`): `table` is `'<table>'` in the writer's
 * schema or `'<schema>.<table>'`; an insert may name itself with `ref`, and any op's value
 * `'$<ref>'` — earlier or later in the batch — is replaced by that row's id. */
export interface DomainTransactionOpLike {
  op: 'insert' | 'update' | 'delete';
  table: string;
  ref?: string;
  values?: Record<string, unknown>;
  id?: string;
  expectedVersion?: number;
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
  access(): Promise<AccessData>;
  query(spec: DomainQueryLike): Promise<Record<string, unknown>[]>;
  /** The total under the same `filter`, `search` and `deleted` mode a query takes. */
  count(filter?: DomainQueryLike['filter'], search?: string, deleted?: DomainDeletedMode): Promise<number>;
  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]>;
  /** The rows as a frame with the writer attached — the one collection a source holds. */
  frame(spec: DomainQueryLike): Promise<DomainFrameLike>;
  /** A row's history, oldest first — the platform's `client.audit`. */
  audit?(id: string): Promise<AuditEntryLike[]>;
  /** Brings a soft-deleted row back (`~is_deleted` off, the `'undelete'` audit op); refused with
   * `validation` where the row refers to a deleted parent. A backend that does not declare it
   * cannot answer a `deleted` query either — `DomainSource` refuses one over it. */
  restore?(id: string): Promise<void>;
}

export interface DomainBackend {
  table(address: string): Promise<DomainTableLike>;
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
