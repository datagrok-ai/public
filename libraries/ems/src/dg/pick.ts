/* `domains.pick`: u2's `PickInput` over a domain table, querying it by its name column through the
   domain seam, so it works over the memory backend as it does over the server. */
import {emsBackends} from '../core/domain-backend.js';
import {InputOptions, text, TypeAhead, Filters, Rows, PickInput} from '@datagrok-libraries/u2';
import type {ObjectRenderer, DomainCondition, DomainConditionNode, DomainConditionTree,
  FilterSchema} from '@datagrok-libraries/u2';
import type {DomainQueryLike, DomainTableInfoLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';

export {PickInput};
export type {PickInputOptions} from '@datagrok-libraries/u2';

export interface PickItem {
  id: string;
  name: string;
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
  /** Row presentation; the name as text by default. */
  renderer?: ObjectRenderer<PickItem>;
  /** How many candidates one search shows (default 10). */
  limit?: number;
  debounceMs?: number;
}

export class DomainPick extends PickInput<PickItem> {
  /** The target's plural name once a search has seen its table — the empty row's wording. */
  private static readonly _plural = new Map<string, string>();

  constructor(readonly table: string, options: DomainPickOptions = {}) {
    const {filter, params, siblings, placeholder, renderer, limit, debounceMs, ...rest} = options;
    const listItem = renderer?.listItem?.bind(renderer);
    let unbound: string[] = [];
    const caption = (name: string) => siblings?.properties.find((p) => p.name === name)?.friendlyName ?? name;
    super({
      ...rest,
      typeAhead: () => new TypeAhead<PickItem>({
        source: (query) => {
          const bound = DomainPick.bind(filter, params?.());
          unbound = bound.unbound;
          return unbound.length > 0 ? Promise.resolve([]) :
            DomainPick.search(table, query, {filter: bound.filter, limit});
        },
        itemText: (item) => item.name,
        render: listItem,
        placeholder,
        debounceMs,
        openOnFocus: true,
        // an empty table is not a dead end: say what there is none of yet
        emptyText: (query) => unbound.length > 0 ? `Pick a ${caption(unbound[0]).toLowerCase()} first` :
          query.trim() === '' ?
            `No ${(DomainPick._plural.get(table) ?? table.split('.').pop()!).toLowerCase()} yet` : 'No matches',
      }),
      idOf: (item) => item.id,
      resolve: (id) => DomainPick.resolve(table, id),
    });
    this.root.dataset.u2 = 'domain-pick';
  }

  /** The candidates for `query`: the rows whose name column contains it, in the table's own order
   * (as seeded — statuses, priorities; a declared order is a phase-4 `views` item). */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, t.info.pluralName);
    const nameColumn = DomainPick.nameColumn(t.info);
    const rows = await t.query({
      filter: DomainPick.filter(nameColumn, query, options.filter),
      limit: options.limit ?? 10,
    });
    return rows.map((row) => DomainPick.item(row, nameColumn));
  }

  /** The item behind an id; a draft id is the draft's caption from the live source holding it,
   * no query issued; an id the table does not answer is shown as itself. A backend that declares
   * `resolveNames` answers it there — over the platform that is the registry's resolver, batched
   * and cached, so a formful of refs is one call; every other backend reads the row. */
  static async resolve(table: string, id: string): Promise<PickItem> {
    const draft = Rows.isDraft(id) ? DomainSource.draftOf(id) : undefined;
    if (draft !== undefined) {
      const info = draft.source.schema.info;
      const name = text(draft.row[DomainPick.nameColumn(info)]);
      return {id, name: name === '' ? `New ${info.singularName.toLowerCase() || 'row'}` : name};
    }
    const backend = emsBackends.domain;
    if (backend?.resolveNames !== undefined) {
      const name = (await backend.resolveNames(table, [id]))[id];
      return {id, name: name || id};
    }
    const t = await DomainPick.table(table);
    const [row] = await t.query({filter: {property: 'id', operator: '=', value: id}, limit: 1});
    return row === undefined ? {id, name: id} : DomainPick.item(row, DomainPick.nameColumn(t.info));
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

  static filter(nameColumn: string, query: string,
    extra?: DomainPickOptions['filter']): DomainQueryLike['filter'] {
    const nodes: DomainConditionNode[] = [];
    const q = query.trim();
    if (q !== '')
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

  static item(row: Record<string, unknown>, nameColumn: string): PickItem {
    const name = text(row[nameColumn]);
    return {id: String(row.id), name: name === '' ? String(row.id) : name};
  }

  private static table(address: string): Promise<DomainTableLike> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      return Promise.reject(new Error('no platform backend for domain tables'));
    return backend.table(address);
  }
}
