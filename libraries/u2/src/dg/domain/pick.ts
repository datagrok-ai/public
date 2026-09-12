/* Reference pickers with an `Input<string | null>` face — the value is the id a ref column
   holds, the box shows the name. `PickInput` is the face over any `TypeAhead`; `domainPick` is
   the one over a domain table, querying it by its name column through the domain seam, so it
   works over the memory backend as it does over the server. */
import {Input, InputOptions, labelText} from '../../core/input-base.js';
import {div} from '../../core/elements.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import {text} from '../../core/text.js';
import {TypeAhead} from '../../components/inputs/typeahead.js';
import {Filters} from '../../core/filter/index.js';
import type {DomainCondition, DomainConditionNode, DomainConditionTree} from '../../core/filter/model.js';
import {backends} from '../../sources/backends.js';
import type {DomainQueryLike, DomainTableInfoLike, DomainTableLike} from '../../sources/domain-backend.js';

export interface PickItem {
  id: string;
  name: string;
}

export interface PickInputOptions<T> extends InputOptions<string | null> {
  /** Builds the type-ahead the input is a face of. */
  typeAhead: () => TypeAhead<T>;
  idOf: (item: T) => string;
  /** The item behind a value written from outside — the id a form loaded; null when unknown. */
  resolve: (id: string) => Promise<T | null>;
}

export class PickInput<T> extends Input<string | null, PickInputOptions<T>> {
  private _typeAhead!: TypeAhead<T>;

  constructor(options: PickInputOptions<T>) {
    super(options, null);
    this.root.dataset.u2 = 'pick-input';
  }

  get typeAhead(): TypeAhead<T> {
    return this._typeAhead;
  }

  // runs from the base constructor, so the coupling state lives here rather than in fields the
  // subclass initializers would reset afterwards
  protected createEditor(): HTMLElement {
    const {typeAhead, idOf, resolve} = this.options;
    const picker = typeAhead();
    this._typeAhead = picker;
    // the last pick is held through typing, so blur can put its name back
    let held: T | null = null;
    let resolving: string | null = null;
    let gen = 0;
    const lookUp = async (id: string) => {
      const mine = ++gen;
      resolving = id;
      try {
        const item = await resolve(id);
        if (mine === gen && this.value.peek() === id && item !== null)
          picker.selected.value = item;
      } catch {
        // an unknown id keeps the value and shows nothing
      } finally {
        if (mine === gen)
          resolving = null;
      }
    };
    // a form seeds a string field with '' for an empty cell: empty either way, and never rewritten
    // to null on the input's own account — that write would dirty the row on render
    const empty = (id: string | null) => id === null || id === '';
    this.effect(() => {
      const id = this.value.value;
      const current = picker.selected.peek();
      if ((current === null ? null : idOf(current)) === id)
        return;
      if (empty(id))
        picker.selected.value = null;
      else
        void lookUp(id!);
    });
    this.effect(() => {
      const item = picker.selected.value;
      // a value still awaiting its item is not a clear
      if (item === null && resolving !== null)
        return;
      const id = item === null ? null : idOf(item);
      const current = this.value.peek();
      // a selection dropped by typing keeps the value; one dropped by clearing the text does not
      if (id !== null)
        held = item;
      else if (picker.text.peek() !== '')
        return;
      else
        held = null;
      if (id === null ? !empty(current) : id !== current)
        this.value.value = id;
    });
    const input = picker.root.querySelector('input')!;
    // text that is not a pick never survives a blur: the held pick's name comes back, else the box
    // is cleared — text that looks chosen is worse than none (an Enter still awaiting its
    // candidates is left to land)
    const onBlur = () => {
      if (picker.isPickPending || picker.selected.peek() !== null || picker.text.peek() === '')
        return;
      if (held !== null)
        picker.selected.value = held;
      else
        picker.resync();
    };
    input.addEventListener('blur', onBlur);
    this.own(() => input.removeEventListener('blur', onBlur));
    const clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'u2-input-clear';
    clear.textContent = '✕';
    clear.setAttribute('aria-label', 'Clear');
    const onClear = () => {
      held = null;
      picker.selected.value = null;
      picker.resync();
      this.value.value = null;
      input.focus();
    };
    clear.addEventListener('click', onClear);
    this.own(() => clear.removeEventListener('click', onClear));
    this.effect(() => clear.hidden = empty(this.value.value) && picker.text.value === '');
    return div([picker.root, clear], 'u2-pick-box');
  }
}

export interface DomainPickOptions extends InputOptions<string | null> {
  /** Narrows the candidates: a smart-filter string or a condition tree, AND-ed with the search. */
  filter?: string | DomainCondition | DomainConditionTree;
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
    const {filter, placeholder, renderer, limit, debounceMs, ...rest} = options;
    const listItem = renderer?.listItem?.bind(renderer);
    super({
      ...rest,
      typeAhead: () => new TypeAhead<PickItem>({
        source: (query) => DomainPick.search(table, query, {filter, limit}),
        itemText: (item) => item.name,
        render: listItem,
        placeholder: placeholder ?? `${labelText(options.label) ?? table.split('.').pop()}…`,
        debounceMs,
        openOnFocus: true,
        // an empty table is not a dead end: say what there is none of yet
        emptyText: (query) => query.trim() === '' ?
          `No ${(DomainPick._plural.get(table) ?? table.split('.').pop()!).toLowerCase()} yet` : 'No matches',
      }),
      idOf: (item) => item.id,
      resolve: (id) => DomainPick.resolve(table, id),
    });
    this.root.dataset.u2 = 'domain-pick';
  }

  /** The candidates for `query`: the rows whose name column contains it, by name. */
  static async search(table: string, query: string,
    options: {filter?: DomainPickOptions['filter'], limit?: number} = {}): Promise<PickItem[]> {
    const t = await DomainPick.table(table);
    DomainPick._plural.set(table, t.info.pluralName);
    const nameColumn = DomainPick.nameColumn(t.info);
    const rows = await t.query({
      filter: DomainPick.filter(nameColumn, query, options.filter),
      sort: nameColumn === 'id' ? undefined : nameColumn,
      limit: options.limit ?? 10,
    });
    return rows.map((row) => DomainPick.item(row, nameColumn));
  }

  /** The item behind an id; an id the table does not answer is shown as itself. */
  static async resolve(table: string, id: string): Promise<PickItem> {
    const t = await DomainPick.table(table);
    const [row] = await t.query({filter: {property: 'id', operator: '=', value: id}, limit: 1});
    return row === undefined ? {id, name: id} : DomainPick.item(row, DomainPick.nameColumn(t.info));
  }

  /** The declared name column, else the first business-key column, else the id. */
  static nameColumn(info: DomainTableInfoLike): string {
    return info.nameColumn ?? info.businessKey[0] ?? 'id';
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
    const backend = backends.domain;
    if (backend === undefined)
      return Promise.reject(new Error('no platform backend for domain tables'));
    return backend.table(address);
  }
}

export function domainPick(table: string, options?: DomainPickOptions): DomainPick {
  return new DomainPick(table, options);
}
