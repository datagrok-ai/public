/* Fills the domain seam (`backends.domain`) at dg import, next to the other platform backends:
   one table handle per address over `grok.dapi.domains`, the registry's metadata as plain
   property records, and `frame()` — the path that makes `DomainSource` a frame host: `queryDf`
   for the rows, the js-api `DomainFrameEditor` attached over them as the single writer, the next
   page appended into the same frame. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type {IProperty} from '../../core/property-like.js';
import {backends} from '../../sources/backends.js';
import type {DataFrameLike} from '../../sources/df-bindings.js';
import type {AuditEntryLike, DomainBackend, DomainBatchOptionsLike, DomainBatchReportLike,
  DomainBatchValidationLike, DomainDeletedMode,
  DomainFrameLike, DomainProbeLike, DomainQueryLike, DomainReadScope, DomainSupportLike, DomainTableInfoLike,
  DomainTableLike, DomainTransactionOpLike, DomainTransactionResultLike} from '../../sources/domain-backend.js';
import {DomainBackendError} from '../../sources/domain-backend.js';
import type {EditState} from '../../sources/edit-state.js';
import {EditorEditState} from './editor-state.js';

/** The metadata fields a form or a filter reads off a property, copied one by one: a `DG.Property`
 * is all prototype getters, and `DomainSource` spreads what the backend hands it. */
const COPIED = ['type', 'propertyType', 'semType', 'nullable', 'choices', 'min', 'max', 'step', 'friendlyName',
  'description', 'inputType', 'editor', 'defaultValue', 'format', 'units', 'showSlider', 'showPlusMinus',
  'category'] as const;

/** js-api `DOMAIN_SYSTEM_COLUMNS`, typed as the frame carries them, with the captions a form shows. */
export const SYSTEM_COLUMNS: readonly (readonly [name: string, type: string, caption: string])[] = [
  ['id', 'string', 'Id'], ['version', 'int', 'Version'], ['created_on', 'datetime', 'Created'],
  ['updated_on', 'datetime', 'Updated'], ['author_id', 'string', 'Author']];

/** The phase-3 calls, through a cast: a package compiles these sources against the PUBLISHED
 * `datagrok-api` types (1.27.11 in every plugin's node_modules), so `client.restore` — and
 * `deleted` in `count`'s options, which `dapi.ts` does not take at all yet — would break every
 * plugin build. Typed here exactly as the server answers them, so the cast hides the version and
 * nothing else. The shim STAYS until the js-api release carrying this surface is out; phase 3-7
 * deletes it. */
interface PhaseThreeClient {
  restore(id: string): Promise<{id: string, restored: boolean, version: number}>;
  count(filter: DomainQueryLike['filter'], options: {search?: string, deleted?: DomainDeletedMode}): Promise<number>;
  updateWhere(filter: DomainQueryLike['filter'], values: Record<string, unknown>,
    options?: {limit?: number}): Promise<{updated: number, hasMore: boolean}>;
  pathTo(id: string): Promise<{id: string, name: string}[]>;
  aggregate(spec: {measures: {fn: string, column?: string, as?: string}[]} & DomainReadScope):
    Promise<Record<string, unknown>[]>;
  /** The table's change token: one bump per write transaction that touched its rows. */
  version(): Promise<{seq: number, at: string | null}>;
  /** The dry run: `POST …/{table}/batch?validateOnly=true`, which answers a verdict per row and
   * writes nothing. */
  batch(rows: Record<string, unknown>[], options: DomainBatchOptionsLike & {validateOnly: true}):
    Promise<DomainBatchValidationLike>;
  access(): Promise<{can: Record<string, boolean>, fields: Record<string, string>,
    support?: DomainSupportLike}>;
}

/** {@link DgDomainBackend.saveAll}'s share of the same debt: `DomainSession.lastRefusal` is not
 * in the published types either. */
interface PhaseThreeSession {
  lastRefusal: {message: string, editor: DG.DomainFrameEditor | null, opIndex?: number} | null;
}

export class DgDomainBackend implements DomainBackend {
  private readonly _tables = new Map<string, Promise<DgDomainTable>>();

  table(address: string): Promise<DomainTableLike> {
    let loading = this._tables.get(address);
    if (loading === undefined) {
      loading = DgDomainTable.load(address).catch((e) => {
        this._tables.delete(address);
        throw e;
      });
      this._tables.set(address, loading);
    }
    return loading;
  }

  resolveNames(table: string, ids: readonly string[]): Promise<Record<string, string | null>> {
    return grok.dapi.domains.registry.resolveNames(table, [...ids]);
  }

  /** Drops every handle, so the next source re-reads the registry — after a grant change. */
  invalidate(): void {
    this._tables.clear();
  }

  /** Every writer's batch as ONE `/transaction` through the js-api `DomainSession`, which owns
   * the conflict dialog, the validation mapping and the retry cap; quiet — the u2 session says
   * what was saved. */
  async saveAll(edits: EditState[]): Promise<boolean> {
    const states = edits as EditorEditState[];
    for (const state of states)
      state.problem = null;
    const session = new DG.DomainSession(states.map((e) => e.editor), {quiet: true});
    try {
      if (await session.save())
        return true;
      // the platform session works out ONE sentence for a refusal ("Ethanol still has 3
      // containers") and balloons it; the u2 status line says it on the state whose editor it
      // names, and on every one for a refusal the batch as a whole answers for
      const refusal = (session as unknown as PhaseThreeSession).lastRefusal;
      for (const state of states) {
        if (refusal !== null && (refusal.editor === null || refusal.editor === state.editor))
          state.problem = refusal.message;
      }
      return false;
    } finally {
      session.dispose();
    }
  }
}

export class DgDomainTable implements DomainTableLike {
  readonly properties: IProperty[];
  /** The six optional members, installed only where {@link support} says the table has them:
   * a control checks `=== undefined` and refuses by name, and never has to guess. */
  readonly audit?: DomainTableLike['audit'];
  readonly restore?: DomainTableLike['restore'];
  readonly updateWhere?: DomainTableLike['updateWhere'];
  readonly ancestors?: DomainTableLike['ancestors'];
  readonly probe?: DomainTableLike['probe'];
  readonly batch?: DomainTableLike['batch'];
  readonly validate?: DomainTableLike['validate'];

  private constructor(readonly address: string, readonly client: DG.DomainTableClient,
    readonly registryProperties: DG.Property[], readonly tableInfo: DG.DomainTableInfo,
    readonly support: DomainSupportLike) {
    this.properties = DgDomainTable.withSystem(registryProperties.map((p) => DgDomainTable.property(p)),
      support.systemColumns);
    if (support.audit)
      this.audit = (id) => this._audit(id);
    if (support.restore)
      this.restore = (id) => this._restore(id);
    if (support.writes) {
      this.updateWhere = (filter, values, options) => this._phase3.updateWhere(filter, values, options);
      this.batch = (rows, options) => this.client.batch(rows, options as DG.DomainBatchOptions) as
        Promise<DomainBatchReportLike>;
      // the dry run is the same endpoint, and it is there wherever the commit is
      this.validate = (rows, options) => this._phase3.batch(rows, {...options, validateOnly: true});
    }
    if (support.ancestors)
      this.ancestors = (id) => this._phase3.pathTo(id);
    if (support.probe)
      this.probe = (scope) => this._probe(scope);
  }

  /** The handle keeps `support` for its lifetime — one registry generation, exactly as the
   * js-api's access cache is; `can`/`fields` are read fresh at every {@link frame}. A server that
   * does not declare it is one boundary failure, never a permissive default: guessing what a
   * table can do is what the declaration replaces. */
  static async load(address: string): Promise<DgDomainTable> {
    const client = grok.dapi.domains.table(address);
    const registry = grok.dapi.domains.registry;
    const [properties, info, access] = await Promise.all([registry.rowProperties(address),
      registry.tableInfo(address), (client as unknown as PhaseThreeClient).access()]);
    if (access.support === undefined) {
      throw new DomainBackendError('unsupported',
        `${address}: the server did not declare its support — restart Datlas on this branch`);
    }
    return new DgDomainTable(address, client, properties, info, access.support);
  }

  /** The system columns every row carries, ahead of the declared ones — `names` is what this
   * table physically has (`support.systemColumns`), so a registration declaring a subset lists
   * only that; the registry describes declared columns only. */
  static withSystem(declared: IProperty[], names: readonly string[]): IProperty[] {
    const own = new Set(declared.map((p) => p.name));
    const system = SYSTEM_COLUMNS.filter(([name]) => !own.has(name) && names.includes(name))
      .map(([name, type, caption]): IProperty => ({
        name, type, propertyType: type, friendlyName: caption, ...(name === 'author_id' ? {semType: 'User'} : {}),
        get: (row: Record<string, unknown>) => row[name],
      } as IProperty));
    return [...system, ...declared];
  }

  get info(): DomainTableInfoLike {
    return this.tableInfo;
  }

  access(): Promise<DG.DomainAccess> {
    return this.client.access();
  }

  query(spec: DomainQueryLike): Promise<Record<string, unknown>[]> {
    return this.client.query(DgDomainTable.spec(spec));
  }

  count(scope: DomainReadScope = {}): Promise<number> {
    return this._phase3.count(scope.filter, {search: scope.search, deleted: scope.deleted});
  }

  private async _restore(id: string): Promise<void> {
    await this._phase3.restore(id);
  }

  /** ONE aggregate for the poll: how many rows match, and when the newest of them was written —
   * or, over the whole live table, ONE read of the change token, which moves once per write
   * transaction and costs no scan at all. `count: -1` says nothing was counted; a source compares
   * the pair against its own previous poll and re-baselines whenever the scope changes. */
  private async _probe(spec: DomainReadScope = {}): Promise<DomainProbeLike> {
    if (DgDomainTable._unscoped(spec)) {
      const version = await this._phase3.version();
      return {count: -1, last: String(version.seq)};
    }
    const rows = await this._phase3.aggregate({
      measures: [{fn: 'count'}, {fn: 'max', column: 'updated_on', as: 'last'}],
      ...(spec.filter === undefined ? {} : {filter: spec.filter}),
      ...(spec.search === undefined ? {} : {search: spec.search}),
      ...(spec.deleted === undefined ? {} : {deleted: spec.deleted}),
    });
    const row = rows[0] ?? {};
    return {count: Number(row.count ?? 0), last: row.last == null ? null : String(row.last)};
  }

  /** Whether a read selects nothing at all — the whole live table, which the change token answers
   * for. An empty condition tree is what an empty filter builder compiles to. */
  private static _unscoped(scope: DomainReadScope): boolean {
    const filter = scope.filter;
    return (filter === undefined || filter === '' || (Array.isArray(filter) && filter.length === 0)) &&
      (scope.search === undefined || scope.search === '') &&
      (scope.deleted === undefined || scope.deleted === 'exclude');
  }

  private get _phase3(): PhaseThreeClient {
    return this.client as unknown as PhaseThreeClient;
  }

  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    return grok.dapi.domains.transaction(this.client.schema, ops as DG.DomainTransactionOp[]);
  }

  /** The row's history as the seam shapes it: `id` is the ROW (the memory backend's key), the
   * platform's audit sequence number is not kept. */
  private async _audit(id: string): Promise<AuditEntryLike[]> {
    const entries = await this.client.audit(id);
    return entries.map((e) => ({...e, id, tx_id: e.tx_id === null ? '' : String(e.tx_id)}));
  }

  /** The frame with the editor attached over it; `queryDf` already stamps the `~can_*` columns out
   * of every export, and an appended page reuses the frame's columns. */
  async frame(spec: DomainQueryLike): Promise<DomainFrameLike> {
    const query = DgDomainTable.spec(spec);
    const [df, access] = await Promise.all([this.client.queryDf(query), this.client.access()]);
    // quiet: the u2 buttons say what was saved; the editor keeps its error and conflict dialogs.
    // A frame over deleted rows is read-only until they are restored — the writer is handed the
    // same upper bound `DomainSource.access` publishes, so it and the controls agree.
    const editor = await DG.DomainFrameEditor.attach(df, this.client,
      {query, access: DgDomainTable.narrowed(access, spec.deleted), quiet: true});
    const edit = new EditorEditState(editor);
    return {
      df: df as unknown as DataFrameLike,
      edit,
      append: async (page) => {
        const rows = await this.client.queryDf(DgDomainTable.spec(page));
        df.append(rows, true);
        return rows.rowCount;
      },
      dispose: () => edit.dispose(),
    };
  }

  /** The seam's query is the js-api spec: same filter tree, same paging keys. */
  static spec(spec: DomainQueryLike): DG.DomainQuerySpec {
    return spec as DG.DomainQuerySpec;
  }

  /** The access a frame's writer is built under: no edit and no insert over deleted rows. */
  static narrowed(access: DG.DomainAccess, deleted: DomainDeletedMode | undefined): DG.DomainAccess {
    return deleted === undefined || deleted === 'exclude' ? access :
      {...access, can: {...access.can, edit: false, insert: false}};
  }

  /** A plain record over a registry property, get/set over a row record — the memory backend's
   * own shape, so a form built from either reads and writes rows the same way. */
  static property(p: DG.Property): IProperty {
    const name = p.name;
    const record: Record<string, unknown> = {
      name,
      get: (row: Record<string, unknown>) => row[name],
      set: (row: Record<string, unknown>, value: unknown) => row[name] = value,
    };
    for (const key of COPIED) {
      const value = p[key];
      if (value !== null && value !== undefined && value !== '')
        record[key] = value;
    }
    return record as IProperty;
  }
}

backends.domain = new DgDomainBackend();
