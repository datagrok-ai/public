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
import type {DomainBackend, DomainFrameLike, DomainQueryLike, DomainTableInfoLike, DomainTableLike,
  DomainTransactionOpLike, DomainTransactionResultLike} from '../../sources/domain-backend.js';
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

  /** Drops every handle, so the next source re-reads the registry — after a grant change. */
  invalidate(): void {
    this._tables.clear();
  }
}

export class DgDomainTable implements DomainTableLike {
  readonly properties: IProperty[];

  private constructor(readonly address: string, readonly client: DG.DomainTableClient,
    readonly registryProperties: DG.Property[], readonly tableInfo: DG.DomainTableInfo) {
    this.properties = DgDomainTable.withSystem(registryProperties.map((p) => DgDomainTable.property(p)));
  }

  static async load(address: string): Promise<DgDomainTable> {
    const client = grok.dapi.domains.table(address);
    const registry = grok.dapi.domains.registry;
    const [properties, info] = await Promise.all([registry.rowProperties(address), registry.tableInfo(address)]);
    return new DgDomainTable(address, client, properties, info);
  }

  /** The system columns every row carries, ahead of the declared ones — as the memory backend
   * lists them, since the registry describes declared columns only. */
  static withSystem(declared: IProperty[]): IProperty[] {
    const own = new Set(declared.map((p) => p.name));
    const system = SYSTEM_COLUMNS.filter(([name]) => !own.has(name)).map(([name, type, caption]): IProperty => ({
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

  count(filter?: DomainQueryLike['filter']): Promise<number> {
    return this.client.count(filter as DG.DomainFilter | undefined);
  }

  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    return grok.dapi.domains.transaction(this.client.schema, ops as DG.DomainTransactionOp[]);
  }

  /** The frame with the editor attached over it; `queryDf` already stamps the `~can_*` columns out
   * of every export, and an appended page reuses the frame's columns. */
  async frame(spec: DomainQueryLike): Promise<DomainFrameLike> {
    const query = DgDomainTable.spec(spec);
    const [df, access] = await Promise.all([this.client.queryDf(query), this.client.access()]);
    // quiet: the u2 buttons say what was saved; the editor keeps its error and conflict dialogs
    const editor = await DG.DomainFrameEditor.attach(df, this.client, {query, access, quiet: true});
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
