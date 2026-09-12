/* The platform seam of `DomainSource` (GOAL "What it looks like"): one table handle answering
   its shape, its access, its rows and its writes. dg fills `backends.domain` over
   `grok.dapi.domains`; `MemoryDomainBackend` answers it from a schema.json in memory, so the
   gallery, the headless tests and an agent iterating on an app need no server. What both must
   agree on is written down in docs/domain-backend-contract.md. */
import type {IProperty} from '../core/property-like.js';
import type {AccessData} from '../core/access.js';
import type {DomainCondition, DomainConditionTree} from '../core/filter/model.js';
import type {DataFrameLike} from './df-bindings.js';
import type {EditState} from './edit-state.js';

export interface DomainTableInfoLike {
  nameColumn: string | null;
  singularName: string;
  pluralName: string;
  businessKey: string[];
}

/** The `DomainQuerySpec` subset a source issues: a smart-filter string or the canonical
 * condition tree, `'col,!col'` ordering, paging, and whether the rows carry the per-row access
 * columns (`Access.ROW_COLUMNS`). */
export interface DomainQueryLike {
  filter?: string | DomainCondition | DomainConditionTree;
  sort?: string;
  columns?: string[];
  limit?: number;
  offset?: number;
  withAccess?: boolean;
}

/** One `/transaction` op (js-api `DomainTransactionOp`): an insert may name itself with `ref`,
 * and a later op's value `'$<ref>'` is replaced by that row's id. */
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
  count(filter?: DomainQueryLike['filter']): Promise<number>;
  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]>;
  /** The rows as a frame with the writer attached — the one collection a source holds. */
  frame(spec: DomainQueryLike): Promise<DomainFrameLike>;
}

export interface DomainBackend {
  table(address: string): Promise<DomainTableLike>;
}

/** What a backend throws: the server's error family by `code` (`'version-conflict'`,
 * `'validation'`, `'not-found'`, `'bad-ref'`), the message for the user. */
export class DomainBackendError extends Error {
  constructor(readonly code: string, message: string) {
    super(message);
    this.name = 'DomainBackendError';
  }
}
