/* `domains.pick`: u2's `PickInput` over a domain table, searching it as the table's own search does
   through the domain seam, so it works over the memory backend as it does over the server. */
import {emsBackends, nounOf} from '../core/domain-backend.js';
import {InputOptions, text, div, span, TypeAhead, Filters, Rows, PickInput} from '@datagrok-libraries/u2';
import type {ObjectRenderer, DomainCondition, DomainConditionNode, DomainConditionTree,
  FilterSchema} from '@datagrok-libraries/u2';
import type {DomainTableInfoLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';
import type {DomainTable} from './facade.js';

export {PickInput};
export type {PickInputOptions} from '@datagrok-libraries/u2';

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

export class DomainPick extends PickInput<PickItem> {
  static host: PickHost | undefined;
  private static readonly _short = new WeakMap<DomainTable, Promise<boolean>>();
  /** The target's plural name once a search has seen its table — the empty row's wording. */
  private static readonly _plural = new Map<string, string>();

  constructor(readonly table: string, options: DomainPickOptions = {}) {
    const {filter, params, siblings, placeholder, renderer, limit, debounceMs, ...rest} = options;
    const listItem = renderer?.listItem?.bind(renderer);
    let unbound: string[] = [];
    let shown: {count: number, total: number} | null = null;
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
            shown = {count: page.items.length, total: page.total};
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
          footer: () => shown === null || shown.total <= shown.count ? null :
            `${shown.count} of ${shown.total} — type to narrow`,
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

  /** Whether the whole table fits one popup of ten; counted once per table and backend. */
  private static async isShort(table: string): Promise<boolean> {
    if (DomainPick.host === undefined)
      return false;
    const handle = await DomainPick.host.handle(table);
    let short = DomainPick._short.get(handle);
    if (short === undefined)
      DomainPick._short.set(handle, short = handle.table.count().then((n) => n <= 10));
    return short;
  }

  /** The candidates for `query`, sorted by what they show. */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    return (await DomainPick.page(table, query, options)).items;
  }

  /** One popup's worth, by the table's own search, sorted by caption — a row whose key is typed
   * (`482`, `DG-000482`) first — and how many rows match in all. */
  private static async page(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<{items: PickItem[], total: number}> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, t.info.pluralName);
    const caption = await DomainPick.host?.caption(table) ?? null;
    const limit = options.limit ?? 10;
    const q = query.trim();
    const extra = DomainPick.filter(undefined, '', options.filter);
    const scope = q !== '' && t.info.searchableColumns.length > 0 ? {search: q, filter: extra} :
      {filter: DomainPick.filter(DomainPick.nameColumn(t.info), q, options.filter)};
    const sort = t.info.nameColumn ?? (t.info.businessKey.length === 1 ? t.info.businessKey[0] : undefined);
    const rows = await t.query({...scope, sort, limit: limit + 1});
    const total = rows.length > limit ? await t.count(scope) : rows.length;
    const items = rows.slice(0, limit).map((row) => DomainPick.item(row, t.info, caption))
      .sort((a, b) => a.name.localeCompare(b.name, undefined, {numeric: true, sensitivity: 'base'}));
    const key = DomainPick.keyCondition(t, q);
    if (key === undefined)
      return {items, total};
    const [row] = await t.query({filter: extra === undefined ? [key] : [key, 'and', ...extra], limit: 1});
    return row === undefined ? {items, total} :
      {items: [DomainPick.item(row, t.info, caption), ...items.filter((item) => item.id !== String(row.id))], total};
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
      return {id, name: name || id};
    }
    const t = await DomainPick.table(table);
    const [row] = await t.query({filter: {property: 'id', operator: '=', value: id}, limit: 1});
    return row === undefined ? {id, name: id} : DomainPick.item(row, t.info, caption);
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

  /** A query that is the table's single numeric key, bare or in its format (`DG-000482`), as a condition. */
  private static keyCondition(t: DomainTableLike, query: string): DomainCondition | undefined {
    const key = t.info.businessKey.length === 1 ? t.info.businessKey[0] : undefined;
    const prop = t.properties.find((p) => p.name === key);
    const type = prop?.propertyType ?? prop?.type;
    if (prop === undefined || (type !== 'int' && type !== 'bigint') || query === '')
      return undefined;
    const prefix = /^[^0#]*/.exec(prop.format ?? '')![0];
    const digits = prefix !== '' && query.toUpperCase().startsWith(prefix.toUpperCase()) ?
      query.slice(prefix.length) : query;
    return /^\d{1,15}$/.test(digits) ? {property: key!, operator: '=', value: Number(digits)} : undefined;
  }

  private static table(address: string): Promise<DomainTableLike> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      return Promise.reject(new Error('no platform backend for domain tables'));
    return backend.table(address);
  }
}
