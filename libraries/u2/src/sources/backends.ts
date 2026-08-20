/* The platform seam of the v1 data sources (Q5): dg fills it at import, headless tests fake it.
   Type skeleton only — the source implementations and the dg filling arrive with WO-14. */
import type {PropertyLike} from '../core/property-like.js';
import type {ObservableLike} from '../core/widget-like.js';
import type {DapiPagerSourceLike} from '../dg/entities/dapi-source.js';
import type {DataFrameLike} from './df-bindings.js';

export interface FuncDescriptorLike {
  name: string;
  /** 'table-query' — a visual query — is the DD9 read-only class: the only kind that is known to
   * read, and the only one the row cap can reach. 'query' is a hand-written one, which the
   * connector executes as written, DML and DDL included. */
  kind: 'table-query' | 'query' | 'script' | 'function' | 'unknown';
  inputs: PropertyLike[];
  outputs: PropertyLike[];
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

export interface SourceBackends {
  funcRunner?: FuncRunnerLike;
  /** Fresh source per call — HttpDataSource's fluent setters mutate (the dapiPager rule). */
  dapi?: (collection: string) => DapiPagerSourceLike<unknown> | null;
  dapiFind?: (collection: string, id: string) => Promise<unknown | null>;
  workspace?: WorkspaceLike;
  /** Sample rows → a DataFrameLike (DG.DataFrame.fromObjects in the platform). */
  tableFromRows?: (rows: object[]) => DataFrameLike;
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
