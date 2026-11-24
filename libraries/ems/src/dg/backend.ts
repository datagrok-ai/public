import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Access, Filters, Rows} from '@datagrok-libraries/u2';
import type {IProperty, DataFrameLike} from '@datagrok-libraries/u2';
import type {AuditEntryLike, DomainBackend, DomainBatchReportLike, DomainDeletedMode, DomainFrameLike, DomainProbeLike,
  DomainQueryLike, DomainReadScope, DomainSupportLike, DomainTableInfoLike, DomainTableLike, DomainTransactionOpLike,
  DomainTransactionResultLike} from '../core/domain-backend.js';
import {emsBackends, DomainBackendError, isUnscoped, nounOf, readableInfo} from '../core/domain-backend.js';
import type {EditState} from '../core/edit-state.js';
import {EditorEditState} from './editor-state.js';

/** The metadata fields a form or a filter reads off a property, copied one by one: a `DG.Property`
 * is all prototype getters, and `DomainSource` spreads what the backend hands it. */
const COPIED = ['type', 'propertyType', 'semType', 'nullable', 'choices', 'min', 'max', 'step', 'friendlyName',
  'description', 'inputType', 'editor', 'defaultValue', 'format', 'units', 'showSlider', 'showPlusMinus',
  'category'] as const;

/** js-api `DOMAIN_SYSTEM_COLUMNS` as the frame types them, with a form's captions; `satisfies`
 * fails the build when js-api's list changes. */
const SYSTEM_SHAPES = {
  id: ['string', 'Id'], version: ['int', 'Version'], created_on: ['datetime', 'Created'],
  updated_on: ['datetime', 'Updated'], author_id: ['string', 'Author'],
} satisfies Record<typeof DG.DOMAIN_SYSTEM_COLUMNS[number], readonly [type: string, caption: string]>;

export const SYSTEM_COLUMNS: readonly (readonly [name: string, type: string, caption: string])[] =
  Object.entries(SYSTEM_SHAPES).map(([name, [type, caption]]) => [name, type, caption]);

/** Why `parent` cannot become the parent of the rows being moved: `among` finds one of them in it and its
 * ancestors, once the pending `moved` (row id → new parent) have landed. */
export async function hierarchyCycle(table: DomainTableLike, parent: unknown,
  among: (ids: string[]) => boolean | Promise<boolean>, moved: ReadonlyMap<unknown, unknown> = new Map(),
  name?: unknown): Promise<string | null> {
  const jumped = new Set<unknown>();
  for (let from = parent; ;) {
    const chain = typeof from !== 'string' || from === '' || Rows.isDraft(from) || table.ancestors === undefined ? [] :
      [from, ...(await table.ancestors(from)).map((a) => a.id).reverse()];
    const at = chain.findIndex((id) => moved.has(id));
    const walked = at < 0 ? chain : chain.slice(0, at + 1);
    if (walked.length > 0 && await among(walked)) {
      return `Cannot move ${name == null ? `a ${nounOf(table.info, true)}` : `"${name}"`} ` +
        'under itself or one of its descendants';
    }
    if (at < 0 || jumped.has(chain[at]))
      return null;
    jumped.add(chain[at]);
    from = moved.get(chain[at]);
  }
}

/** Fills `emsBackends.domain` at dg import over `grok.dapi.domains`; its `frame()` makes `DomainSource` a frame
 * host: `queryDf` for the rows, the js-api `DomainFrameEditor` over them as the single writer. */
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

  dropAccess(): void {
    grok.dapi.domains.invalidateUiCaches();
    this.invalidate();
  }

  /** A refusal the server answered with 403: by its code, or by the message prefix where js-api records none. */
  static isForbidden(refusal: {message?: string, code?: string} | null): boolean {
    if (refusal == null)
      return false;
    if (refusal.code !== undefined)
      return refusal.code === 'forbidden';
    const message = refusal.message ?? '';
    return message === 'forbidden' || message.startsWith('Insufficient privileges to ') ||
      message.startsWith('You don\'t have ');
  }

  /** Every writer's batch as ONE `/transaction` through the js-api `DomainSession`, which owns
   * the conflict dialog, the validation mapping and the retry cap; quiet — the u2 session says
   * what was saved. */
  async saveAll(edits: EditState[]): Promise<boolean> {
    const states = edits as EditorEditState[];
    for (const state of states) {
      state.problem = null;
      state.problemCode = null;
      state.retry?.();
    }
    const moves = states.flatMap((state) => {
      const {editor, df, table} = state;
      const column = table?.info.parentColumn;
      const name = table?.info.nameColumn;
      return column == null ? [] : [...Array(df.rowCount).keys()]
        .filter((i) => editor.stateOf(i) === 'modified' && editor.isChanged(i, column))
        .map((i) => ({state, id: df.get('id', i), parent: df.get(column, i),
          name: name == null ? null : df.get(name, i)}));
    });
    const moved = new Map(moves.map((m) => [m.id, m.parent]));
    for (const {state, id, parent, name} of moves) {
      state.problem = await hierarchyCycle(state.table!, parent, (ids) => ids.includes(id as string), moved, name);
      if (state.problem !== null)
        return false;
    }
    const before = states.map((state) => state.errorCells?.());
    const session = new DG.DomainSession(states.map((e) => e.editor), {quiet: true});
    try {
      if (await session.save())
        return true;
      // the platform session works out ONE sentence for a refusal ("Ethanol still has 3
      // containers") and balloons it; the u2 status line says it on the state whose editor it
      // names, and on every one for a refusal the batch as a whole answers for
      const refusal = session.lastRefusal;
      const sent = refusal !== null && refusal.message !== DG.DomainSession.CELL_ERRORS;
      states.forEach((state, i) => state.markRefused?.(sent ? before[i]! : state.errorCells()));
      for (const state of states) {
        if (refusal !== null && (refusal.editor === null || refusal.editor === state.editor))
          state.problem = refusal.message;
        state.problemCode = DgDomainBackend.isForbidden(refusal) ? 'forbidden' : refusal?.code ?? null;
      }
      return false;
    } finally {
      session.dispose();
    }
  }
}

export class DgDomainTable implements DomainTableLike {
  readonly properties: IProperty[];
  /** The seven optional members, installed only where {@link support} says the table has them:
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
    readonly support: DomainSupportLike, private _token: boolean) {
    this.properties = DgDomainTable.withSystem(registryProperties.map((p) => DgDomainTable.property(p)),
      support.systemColumns);
    if (support.audit)
      this.audit = (id, page) => this._audit(id, page);
    if (support.restore)
      this.restore = (id) => this._restore(id);
    if (support.updateWhere) {
      this.updateWhere = (filter, values, options) =>
        this.client.updateWhere(DgDomainTable.bound(filter) as DG.DomainFilter, values, options);
    }
    if (support.writes) {
      this.batch = (rows, options) => this.client.batch(rows, options as DG.DomainBatchOptions) as
        Promise<DomainBatchReportLike>;
    }
    if (support.batch.validate) {
      this.validate = (rows, options) => this.client.batch(rows,
        {...options, validateOnly: true} as DG.DomainBatchOptions & {validateOnly: true});
    }
    if (support.ancestors)
      this.ancestors = (id) => this.client.pathTo(id);
    if (support.probe)
      this.probe = (scope) => this._probe(scope);
  }

  /** The handle keeps `support` for its lifetime — one registry generation, exactly as the
   * js-api's access cache is; `can`/`fields` are read fresh at every {@link frame}. A server that
   * does not declare it is one boundary failure, never a permissive default: guessing what a
   * table can do is what the declaration replaces. One member is checked — `support.batch`, the
   * newest storage switch — and stands for the whole contract generation. */
  static async load(address: string): Promise<DgDomainTable> {
    const client = grok.dapi.domains.table(address);
    const registry = grok.dapi.domains.registry;
    const [properties, info, access] = await Promise.all([registry.rowProperties(address),
      registry.tableInfo(address), client.access()]);
    if (access.support?.batch == null) {
      throw new DomainBackendError('unsupported',
        `${address}: the server declares no storage switches (support.batch) — restart Datlas on this branch`);
    }
    // the change token is read behind the table-level View grant: a reader who reaches rows
    // through row grants only polls the aggregate over what it may see
    const readable = DgDomainTable.readable(access);
    return new DgDomainTable(address, client, properties.filter((p) => readable(p.name)),
      readableInfo(info, properties, readable), access.support,
      access.support.version && access.can.view);
  }

  /** Whether the caller may read a column: only `access().fields` leaves a restricted one out. */
  static readable(access: DG.DomainAccess): (name: string) => boolean {
    const fields = Access.from(access);
    return (name) => fields.field(name) !== 'hidden';
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
    return this.client.count({...scope, filter: DgDomainTable.bound(scope.filter)} as DG.DomainReadScope);
  }

  private async _restore(id: string): Promise<void> {
    await this.client.restore(id);
  }

  // over the whole live table the change token answers with no scan; `count: -1` says nothing was counted
  private async _probe(spec: DomainReadScope = {}): Promise<DomainProbeLike> {
    if (this._token && isUnscoped(spec)) {
      try {
        const version = await this.client.version();
        return {count: -1, last: String(version.seq)};
      } catch {
        // the grant behind the token went away mid-session: the aggregate answers from here on,
        // instead of the failure counter stopping the poll
        this._token = false;
      }
    }
    const rows = await this.client.aggregate({
      measures: [{fn: 'count'}, {fn: 'max', column: 'updated_on', as: 'last'}],
      ...(spec.filter === undefined ? {} : {filter: DgDomainTable.bound(spec.filter) as DG.DomainFilter}),
      ...(spec.search === undefined ? {} : {search: spec.search}),
      ...(spec.deleted === undefined ? {} : {deleted: spec.deleted}),
    });
    const row: Record<string, unknown> = rows[0] ?? {};
    return {count: Number(row.count ?? 0), last: row.last == null ? null : String(row.last)};
  }

  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    return grok.dapi.domains.transaction(this.client.schema, ops as DG.DomainTransactionOp[]);
  }

  // the seam's `id` is the row (the memory backend's key); the platform's audit sequence number is `seq`
  private async _audit(id: string, page?: {limit?: number, before?: number}): Promise<AuditEntryLike[]> {
    const entries = await this.client.audit(id, page);
    return entries.map((e) => ({...e, id, seq: e.id, tx_id: e.tx_id === null ? '' : String(e.tx_id)}));
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
      {query, access: DgDomainTable.narrowed(access, spec.deleted, this.support.transaction), quiet: true});
    const edit = new EditorEditState(editor, this);
    return {
      df: df as unknown as DataFrameLike,
      edit,
      append: async (page) => {
        const rows = await this.client.queryDf(DgDomainTable.spec(page));
        const read = rows.rowCount;
        const ids = df.col('id');
        const theirs = rows.col('id');
        if (ids !== null && theirs !== null) {
          const held = new Set(ids.toList());
          rows.rows.removeWhereIdx((i) => held.has(theirs.get(i)));
        }
        df.append(rows, true);
        return read;
      },
      dispose: () => edit.dispose(),
    };
  }

  /** The seam's query is the js-api spec: same filter tree, same paging keys. */
  static spec(spec: DomainQueryLike): DG.DomainQuerySpec {
    return {...spec, filter: DgDomainTable.bound(spec.filter)} as DG.DomainQuerySpec;
  }

  /** `$me` bound to the caller: the query an app keeps, its address and its recent searches say `$me`. */
  static bound(filter: DomainQueryLike['filter']): DomainQueryLike['filter'] {
    const me = {me: grok.shell.user.id};
    if (typeof filter === 'string') {
      const parsed = /\$me\b/.test(filter) ? Filters.parse(filter) : undefined;
      return parsed === undefined || parsed.problems.length > 0 ? filter :
        Filters.format(Filters.bind(parsed.root, me));
    }
    return filter === undefined || !JSON.stringify(filter).includes('"$param":"me"') ? filter :
      Filters.toDomainTree(Filters.bind(Filters.fromDomainTree(filter), me));
  }

  /** The access a frame's writer is built under: no edit and no insert over deleted rows, nothing
   * at all over a table `transaction` does not land on — the bound `DomainSource.access` publishes.
   * The latter closes every field too: the editor lets a row's own `~can_edit` lift the table-level
   * `edit`, and no row grant may reopen what the source cannot save. */
  static narrowed(access: DG.DomainAccess, deleted: DomainDeletedMode | undefined,
    transaction = true): DG.DomainAccess {
    if (!transaction) {
      return {...access, can: {...access.can, edit: false, insert: false, delete: false},
        fields: Object.fromEntries(Object.keys(access.fields).map((c) => [c, 'readonly' as const]))};
    }
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

emsBackends.domain = new DgDomainBackend();
