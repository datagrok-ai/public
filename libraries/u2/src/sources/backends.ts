/* The platform seam (Q5) of the data sources and of the platform defaults the components fall back to — a
   principal's picture, the formatting: dg fills it at import, headless tests fake it. */
import type {IProperty} from '../core/property-like.js';
import type {NamedProperty, ObservableLike} from '../core/widget-like.js';
import type {DapiPagerSourceLike} from '../dg/entities/dapi-source.js';
import type {DataFrameLike} from './df-bindings.js';

export interface FuncDescriptorLike {
  name: string;
  /** 'table-query' — a visual query — is the DD9 read-only class: the only kind that is known to
   * read, and the only one the row cap can reach. 'query' is a hand-written one, which the
   * connector executes as written, DML and DDL included. */
  kind: 'table-query' | 'query' | 'script' | 'function' | 'unknown';
  inputs: NamedProperty[];
  outputs: IProperty[];
}

export interface FuncRunnerLike {
  find(name: string): FuncDescriptorLike | null;
  /** Outputs by declared name; single-output funcs answer {<name>: value}. The platform offers no
   * cancellation for a function call — `FuncCall.cancel()` marks the call and fires an event, but
   * `cancelImpl` is a no-op and nothing on the client sets `isCancelable` — so `abort` reaches
   * whatever the implementation can stop on its own; a superseded run still finishes on the
   * server and the caller drops its answer. */
  run(name: string, params: Record<string, unknown>, abort: AbortSignal,
    options?: {maxRows?: number}): Promise<Record<string, unknown>>;
}

export interface WorkspaceLike {
  table(name: string): DataFrameLike | null;
  tableNames(): string[];
  onTablesChanged: ObservableLike<unknown>;
}

/** The platform's formatting, which `Format` applies. */
export interface FormatBackend {
  /** The date pattern dates read in: the user's Settings › Date format; null or '' for the platform's default. */
  dateFormat(): string | null;
  /** `x` in a platform number format: '#,##0.00', 'compact', 'scientific'. */
  number(x: number, format: string): string;
}

export interface SourceBackends {
  funcRunner?: FuncRunnerLike;
  /** Fresh source per call — HttpDataSource's fluent setters mutate (the dapiPager rule). */
  dapi?: (collection: string) => DapiPagerSourceLike<unknown> | null;
  dapiFind?: (collection: string, id: string) => Promise<unknown | null>;
  workspace?: WorkspaceLike;
  /** Sample rows → a DataFrameLike (DG.DataFrame.fromObjects in the platform). */
  tableFromRows?: (rows: object[]) => DataFrameLike;
  /** The picture before a principal's name in an `AccessGrid` with no `icon` of its own — a principal is a group's
   * id or its name: a user's picture, a group's icon; null where there is none. */
  principalIcon?: (principal: string) => HTMLElement | null;
  format?: FormatBackend;
}

export const backends: SourceBackends = {};

/** The backend a source cannot work without. The scope exists from `super()`, so a source that
 * cannot work disposes itself before throwing rather than leave one behind. */
export function requireBackend<T>(source: {dispose(): void}, backend: T | undefined,
  what: string): T {
  if (backend === undefined) {
    source.dispose();
    throw new Error(`no platform backend for ${what}`);
  }
  return backend;
}
