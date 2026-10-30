import {emsBackends, nounOf} from '../core/domain-backend.js';
import {InputOptions, text, div, span, TypeAhead, Filters, Rows, PickInput} from '@datagrok-libraries/u2';
import type {ObjectRenderer, DomainCondition, DomainConditionNode, DomainConditionTree,
  FilterSchema} from '@datagrok-libraries/u2';
import type {DomainTableInfoLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';
import type {DomainTable} from './facade.js';

export interface PickItem {
  id: string;
  name: string;
  /** A muted second line: the row's name where the caption is an identifier. */
  detail?: string;
}

/** The app's caption for a table's rows (null without an app renderer) and how to open a row; set by
 * `form.ts`, so this module needs no runtime import of the facade. */
export interface PickHost {
  caption(address: string): Promise<((row: Record<string, unknown>) => string) | null>;
  open(address: string, id: string): void;
  handle(address: string): Promise<DomainTable>;
}

export interface DomainPickOptions extends InputOptions<string | null> {
  /** Narrows the candidates: a smart-filter string or a condition tree, AND-ed with the search.
   * A `$name` in the string is a sibling's value, bound through {@link params} at every search. */
  filter?: string | DomainCondition | DomainConditionTree;
  /** The values `$name`s bind to — the row being edited; a `$name` left empty means no
   * candidates and a "Pick a <sibling> first" hint. */
  params?: () => Record<string, unknown>;
  /** The schema the `$name`s belong to — the hint's captions. */
  siblings?: FilterSchema;
  placeholder?: string;
  /** Row presentation; the caption over its detail by default. */
  renderer?: ObjectRenderer<PickItem>;
  /** How many candidates one search shows (default 10). */
  limit?: number;
  debounceMs?: number;
}

/** u2's `PickInput` over a domain table, searching it as the table's own search does through the domain seam, so
 * it works over the memory backend as it does over the server. */
export class DomainPick extends PickInput<PickItem> {
  static host: PickHost | undefined;
  private static readonly _short = new WeakMap<DomainTable, Promise<boolean>>();
  private static readonly _plural = new Map<string, string>();

  constructor(readonly table: string, options: DomainPickOptions = {}) {
    const {filter, params, siblings, placeholder, renderer, limit, debounceMs, ...rest} = options;
    const listItem = renderer?.listItem?.bind(renderer);
    let unbound: string[] = [];
    let more = false;
    const caption = (name: string) => siblings?.properties.find((p) => p.name === name)?.friendlyName ?? name;
    const host = DomainPick.host;
    super({
      ...rest,
      typeAhead: () => {
        const typeAhead: TypeAhead<PickItem> = new TypeAhead<PickItem>({
          source: async (query) => {
            const bound = DomainPick.bind(filter, params?.());
            unbound = bound.unbound;
            if (unbound.length > 0)
              return [];
            const page = await DomainPick.page(table, query, {filter: bound.filter, limit});
            more = page.more;
            // the held value opens the list, so it is never missing from it
            const held = typeAhead.selected.peek();
            return held === null || query.trim() !== '' ? page.items :
              [held, ...page.items.filter((item) => item.id !== held.id)];
          },
          itemText: (item) => item.name,
          keyOf: (item) => item.id,
          render: listItem ?? DomainPick.renderItem,
          placeholder,
          debounceMs,
          openOnFocus: true,
          matchWidth: false,
          footer: () => more ? 'More match — type to narrow' : null,
          // an empty table is not a dead end: say what there is none of yet
          emptyText: (query) => unbound.length > 0 ? `Pick a ${caption(unbound[0]).toLowerCase()} first` :
            query.trim() === '' ?
              `No ${(DomainPick._plural.get(table) ?? table.split('.').pop()!).toLowerCase()} yet` : 'No matches',
        });
        return typeAhead;
      },
      idOf: (item) => item.id,
      resolve: (id) => DomainPick.resolve(table, id),
      open: host === undefined ? undefined : (id) => {
        if (!Rows.isDraft(id))
          host.open(table, id);
      },
    });
    this.root.dataset.u2 = 'domain-pick';
    // a short lookup (a status, a priority) reads as a dropdown
    void DomainPick.isShort(table).then((short) => {
      if (short && !this.scope.isDisposed)
        this.asDropdown();
    }, () => undefined);
  }

  private static async isShort(table: string): Promise<boolean> {
    if (DomainPick.host === undefined)
      return false;
    const handle = await DomainPick.host.handle(table);
    let short = DomainPick._short.get(handle);
    if (short === undefined)
      DomainPick._short.set(handle, short = handle.table.count().then((n) => n <= 10));
    return short;
  }

  /** The candidates for `query`. */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    return (await DomainPick.page(table, query, options)).items;
  }

  private static async page(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<{items: PickItem[], more: boolean}> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, t.info.pluralName);
    const caption = await DomainPick.host?.caption(table) ?? null;
    const limit = options.limit ?? 10;
    const q = query.trim();
    const extra = DomainPick.filter(undefined, '', options.filter);
    const key = DomainPick.keyCondition(t, q);
    const scope = key !== undefined ? {filter: DomainPick.orKey(t.info, q, key, extra)} :
      q !== '' && t.info.searchableColumns.length > 0 ? {search: q, filter: extra} :
      {filter: DomainPick.filter(DomainPick.nameColumn(t.info), q, options.filter)};
    const sort = t.info.nameColumn ?? (t.info.businessKey.length === 1 ? t.info.businessKey[0] : undefined);
    const rows = await t.query({...scope, sort, limit: limit + 1});
    return {items: rows.slice(0, limit).map((row) => DomainPick.item(row, t.info, caption)), more: rows.length > limit};
  }

  /** The item behind an id; a draft id is the draft's caption from the live source holding it,
   * no query issued; an id the table does not answer is shown as itself. A table with an app caption
   * reads the row; otherwise a backend that declares `resolveNames` answers it there — over the
   * platform the registry's resolver, batched and cached; every other backend reads the row. */
  static async resolve(table: string, id: string): Promise<PickItem> {
    const draft = Rows.isDraft(id) ? DomainSource.draftOf(id) : undefined;
    if (draft !== undefined) {
      const info = draft.source.schema.info;
      const name = text(draft.row[DomainPick.nameColumn(info)]);
      return {id, name: name === '' ? `New ${nounOf(info, true)}` : name};
    }
    const caption = await DomainPick.host?.caption(table) ?? null;
    const backend = emsBackends.domain;
    if (caption === null && backend?.resolveNames !== undefined) {
      const name = (await backend.resolveNames(table, [id]))[id];
      // no name from the server is also how an external table answers: its row is read instead
      if (name)
        return {id, name};
    }
    const t = await DomainPick.table(table);
    const row = await DomainPick.rowOf(table, t, id);
    return row === undefined ? {id, name: id} : DomainPick.item(row, t.info, caption);
  }

  private static _pending = new Map<string, Map<string, ((row: Record<string, unknown> | undefined) => void)[]>>();

  /** The row behind an id, read together with every other id of the table asked for in the same tick: a page
   * of links costs one query per table, not one per cell. */
  private static rowOf(address: string, t: DomainTableLike, id: string): Promise<Record<string, unknown> | undefined> {
    let ids = DomainPick._pending.get(address);
    if (ids === undefined) {
      const batch = ids = new Map();
      DomainPick._pending.set(address, batch);
      // the cells waiting on the same lookups resume in one turn: they have all joined by the next microtask
      queueMicrotask(() => {
        DomainPick._pending.delete(address);
        const keys = [...batch.keys()];
        const answer = (rows: Record<string, unknown>[]) => {
          const byId = new Map(rows.map((r) => [String(r.id), r]));
          for (const [key, waiting] of batch) {
            for (const done of waiting)
              done(byId.get(key));
          }
        };
        // `=` over a list is the wire's IN
        t.query({filter: {property: 'id', operator: '=', value: keys}, limit: keys.length})
          .then(answer, () => answer([]));
      });
    }
    return new Promise((done) => ids!.set(id, [...ids!.get(id) ?? [], done]));
  }

  private static renderItem(item: PickItem): HTMLElement {
    const name = span(item.name, 'u2-typeahead-text');
    name.dataset.u2Part = 'item-label';
    return item.detail === undefined ? name :
      div([name, span(item.detail, 'u2-domain-pick-detail')], 'u2-domain-pick-item');
  }

  /** The declared name column, else the first business-key column, else the id. */
  static nameColumn(info: DomainTableInfoLike): string {
    return info.nameColumn ?? info.businessKey[0] ?? 'id';
  }

  /** A string filter with its `$name`s bound to `params` (an empty sibling counts as absent), as a
   * tree; the names still unbound, so the picker can say which sibling to pick first. */
  static bind(extra: DomainPickOptions['filter'], params?: Record<string, unknown>):
    {filter?: DomainCondition | DomainConditionTree, unbound: string[]} {
    if (extra === undefined || (typeof extra === 'string' && extra.trim() === ''))
      return {unbound: []};
    if (typeof extra !== 'string')
      return {filter: extra, unbound: []};
    const values = Object.fromEntries(Object.entries(params ?? {})
      .filter(([, v]) => v !== null && v !== undefined && v !== ''));
    const root = Filters.bind(Filters.parse(extra).root, values);
    const unbound: string[] = [];
    Filters.walk(root, (n) => {
      if (Filters.isGroup(n))
        return;
      for (const v of Array.isArray(n.value) ? n.value : [n.value]) {
        if (Filters.isParam(v) && !unbound.includes(v.param))
          unbound.push(v.param);
      }
    });
    return {filter: unbound.length > 0 ? undefined : Filters.toDomainTree(root), unbound};
  }

  static filter(nameColumn: string | undefined, query: string,
    extra?: DomainPickOptions['filter']): DomainConditionNode[] | undefined {
    const nodes: DomainConditionNode[] = [];
    const q = query.trim();
    if (q !== '' && nameColumn !== undefined)
      nodes.push({property: nameColumn, operator: 'like', value: `%${Filters.escapeLike(q)}%`});
    const more = typeof extra === 'string' ?
      (extra.trim() === '' ? undefined : Filters.toDomainTree(Filters.parse(extra).root)) : extra;
    if (more !== undefined) {
      if (nodes.length > 0)
        nodes.push('and');
      nodes.push(more);
    }
    return nodes.length === 0 ? undefined : nodes;
  }

  /** The app's caption (else the name column, else the id); the name or a searchable value as its
   * detail where the caption is not the name. */
  static item(row: Record<string, unknown>, info: DomainTableInfoLike,
    caption: ((row: Record<string, unknown>) => string) | null = null): PickItem {
    const named = text(row[DomainPick.nameColumn(info)]);
    const name = caption?.(row) || named || String(row.id);
    if (name === named)
      return {id: String(row.id), name};
    const detail = [info.nameColumn, ...info.searchableColumns].filter((c): c is string => c != null)
      .map((c) => text(row[c])).find((v) => v !== '' && v !== name);
    return detail === undefined ? {id: String(row.id), name} : {id: String(row.id), name, detail};
  }

  private static orKey(info: DomainTableInfoLike, q: string, key: DomainConditionNode,
    extra: DomainConditionNode[] | undefined): DomainConditionTree {
    const columns = info.searchableColumns.length > 0 ? info.searchableColumns : [DomainPick.nameColumn(info)];
    const match: DomainConditionNode[] = [key];
    for (const c of columns)
      match.push('or', {property: c, operator: 'like', value: `%${Filters.escapeLike(q)}%`});
    return extra === undefined ? match : [match, 'and', ...extra];
  }

  // `DG-000482` is that number, a part after the prefix a range of numbers, bare digits a substring of the number
  private static keyCondition(t: DomainTableLike, query: string): DomainConditionNode | undefined {
    const column = DomainPick.nameColumn(t.info);
    const prop = t.properties.find((p) => p.name === column);
    const type = prop?.propertyType ?? prop?.type;
    const format = prop?.format ?? '';
    if ((type !== 'int' && type !== 'bigint') || format === '' || query === '')
      return undefined;
    const prefix = /^[^0#]*/.exec(format)![0];
    const width = /^[0#]*/.exec(format.slice(prefix.length))![0].length;
    const anchored = prefix !== '' && query.toUpperCase().startsWith(prefix.toUpperCase());
    const digits = anchored ? query.slice(prefix.length) : query;
    if (!/^\d{1,15}$/.test(digits))
      return undefined;
    if (!anchored) {
      const bare = digits.replace(/^0+/, '');
      return bare === '' ? undefined : {property: column, operator: 'like', value: `%${bare}%`};
    }
    if (digits.length >= width)
      return {property: column, operator: '=', value: Number(digits)};
    const scale = 10 ** (width - digits.length);
    return [{property: column, operator: '>=', value: Number(digits) * scale}, 'and',
      {property: column, operator: '<', value: (Number(digits) + 1) * scale}];
  }

  private static table(address: string): Promise<DomainTableLike> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      return Promise.reject(new Error('no platform backend for domain tables'));
    return backend.table(address);
  }
}
