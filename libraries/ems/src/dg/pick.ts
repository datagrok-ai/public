/* `domains.pick`: u2's `PickInput` over a domain table, searching it as the table's own search does
   through the domain seam, so it works over the memory backend as it does over the server. */
import {emsBackends} from '../core/domain-backend.js';
import {InputOptions, text, div, span, TypeAhead, Filters, Rows, PickInput} from '@datagrok-libraries/u2';
import type {ObjectRenderer, DomainCondition, DomainConditionNode, DomainConditionTree,
  FilterSchema} from '@datagrok-libraries/u2';
import type {DomainTableInfoLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';

export {PickInput};
export type {PickInputOptions} from '@datagrok-libraries/u2';

export interface PickItem {
  id: string;
  name: string;
  /** A muted second line: what the row is called where its caption is an identifier. */
  detail?: string;
}

/** What a picker asks of the layer that knows a table's app: the caption its configured renderer
 * gives a row (null where no app dressed the table — the name column is the caption then), and
 * how to open a row. Installed by the dg layer (`form.ts`); a headless picker has none. */
export interface PickHost {
  caption(address: string): Promise<((row: Record<string, unknown>) => string) | null>;
  open(address: string, id: string): void;
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
  /** The most rows a lookup read as a dropdown holds. */
  static readonly SHORT = 10;
  /** Per backend, so a backend swapped in (a test, a reconnect) counts its own tables. */
  private static readonly _short = new WeakMap<object, Map<string, Promise<boolean>>>();
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
      typeAhead: () => new TypeAhead<PickItem>({
        source: async (query) => {
          const bound = DomainPick.bind(filter, params?.());
          unbound = bound.unbound;
          if (unbound.length > 0)
            return [];
          const page = await DomainPick.page(table, query, {filter: bound.filter, limit});
          shown = {count: page.items.length, total: page.total};
          return page.items;
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
      }),
      idOf: (item) => item.id,
      resolve: (id) => DomainPick.resolve(table, id),
      open: host === undefined ? undefined : (id) => {
        if (!Rows.isDraft(id))
          host.open(table, id);
      },
    });
    this.root.dataset.u2 = 'domain-pick';
    // a short lookup (a status, a priority) reads as a dropdown: a click lists it all, typing has nothing to narrow
    void DomainPick.isShort(table).then((short) => {
      if (!short || this.scope.isDisposed)
        return;
      const input = this.typeAhead.root.querySelector('input')!;
      input.readOnly = true;
      this.root.classList.add('u2-pick-short');
      const onClick = () => {
        if (!this.typeAhead.isOpen.peek())
          this.typeAhead.open();
      };
      input.addEventListener('click', onClick);
      this.own(() => input.removeEventListener('click', onClick));
    }, () => undefined);
  }

  /** Whether the whole table fits one popup — at most {@link SHORT} rows, counted once per table. */
  static isShort(table: string): Promise<boolean> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      return Promise.resolve(false);
    let counted = DomainPick._short.get(backend);
    if (counted === undefined) {
      counted = new Map();
      DomainPick._short.set(backend, counted);
    }
    let short = counted.get(table);
    if (short === undefined) {
      short = DomainPick.table(table).then((t) => t.count()).then((n) => n <= DomainPick.SHORT);
      counted.set(table, short);
    }
    return short;
  }

  /** The candidates for `query`, sorted by what they show. */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    return (await DomainPick.page(table, query, options)).items;
  }

  /** One popup's worth: the rows the table's search finds (its searchable columns, else the name
   * column), the row whose key is written as typed — `482` or `DG-000482` — first, the rest sorted
   * by caption; and how many rows match in all. */
  static async page(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<{items: PickItem[], total: number}> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, t.info.pluralName);
    const caption = await DomainPick.host?.caption(table) ?? null;
    const limit = options.limit ?? 10;
    const q = query.trim();
    const extra = DomainPick.filter(undefined, '', options.filter);
    const scope = q !== '' && t.info.searchableColumns.length > 0 ? {search: q, filter: extra} :
      {filter: DomainPick.filter(DomainPick.nameColumn(t.info), q, options.filter)};
    const rows = await t.query({...scope, sort: DomainPick.sortColumn(t.info), limit: limit + 1});
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
   * no query issued; an id the table does not answer is shown as itself. A table an app dressed
   * reads the row and captions it as the app does; otherwise a backend that declares
   * `resolveNames` answers it there — over the platform that is the registry's resolver, batched
   * and cached, so a formful of refs is one call; every other backend reads the row. */
  static async resolve(table: string, id: string): Promise<PickItem> {
    const draft = Rows.isDraft(id) ? DomainSource.draftOf(id) : undefined;
    if (draft !== undefined) {
      const info = draft.source.schema.info;
      const name = text(draft.row[DomainPick.nameColumn(info)]);
      return {id, name: name === '' ? `New ${info.singularName.toLowerCase() || 'row'}` : name};
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

  /** A candidate row: the caption over a muted line naming the row where the caption does not. */
  static renderItem(item: PickItem): HTMLElement {
    const name = span(item.name, 'u2-typeahead-text');
    name.dataset.u2Part = 'item-label';
    return item.detail === undefined ? name : div([name, span(item.detail, 'u2-domain-pick-detail')], 'u2-domain-pick-item');
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

  /** A row as a candidate: the app's caption (else the name column, else the id), and where that
   * caption is not the row's name, the name — else the first searchable value — as its detail. */
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

  /** The column the server orders candidates by: the name column, else a one-column key. */
  static sortColumn(info: DomainTableInfoLike): string | undefined {
    return info.nameColumn ?? (info.businessKey.length === 1 ? info.businessKey[0] : undefined);
  }

  /** A query that is the table's one numeric key as the table writes it — `482`, or `DG-000482`
   * under the column's `DG-000000` format — as a condition on that key. */
  static keyCondition(t: DomainTableLike, query: string): DomainCondition | undefined {
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
